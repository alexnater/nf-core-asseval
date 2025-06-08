#!/usr/bin/env python

import sys
import re
import math
import logging
from pathlib import Path
from typing import Tuple
from enum import IntEnum
from dataclasses import dataclass
import array as arr
from pysam import TabixFile, VariantFile, FastaFile, asTuple

logger = logging.getLogger(__name__)

BASES = ('A', 'T', 'G', 'C')

class Base(IntEnum):
    A = 0
    T = 1
    G = 2
    C = 3


class Region:
    def __init__(self, chrom: str=None, start: int=None, end: int=None, name: str=None):
        if chrom is None or start is None or end is None:
            self.defined = False
        else:
            self.defined = True
        if self.defined and start > end:
            raise ValueError (f"Start value {start} is larger than end value {end}!")
        self.chrom = chrom
        self.startpos = start
        self.endpos = end
        self.name = name if name else f"{chrom}:{start}-{end}"

    @classmethod
    def from_string(cls, region_str: str):
        tmp = region_str.split(':')
        if len(tmp) != 2:
            raise ValueError (f"Invalid region string {region_str}!")
        positions = tuple(int(x) for x in tmp[1].split('-'))
        if len(positions) != 2:
            raise ValueError (f"Invalid region string {region_str}!")
        return cls(tmp[0], positions[0], positions[1])

    def is_in_interval(self, chrom: str, pos: int) -> bool:
        if not chrom == self.chrom:
            return False
        if pos >= self.startpos and pos <= self.endpos:
            return True
        return False
    
    def __str__(self):
        return f"{self.chrom}_{self.start}-{self.end}"
    
    def __repr__(self):
        return f"{self.chrom}:{self.start}-{self.end}"

    @property    
    def start(self) -> int:
        if self.startpos is None:
            return 1
        else:
            return self.startpos
        
    @property    
    def end(self) -> int:
        return self.endpos

    @property    
    def length(self) -> int:
        if self.defined:
            return self.endpos - self.startpos + 1
        else:
            return None


class Matrix:
    def __init__(self, type: str, nrows: int=0, ncols: int=0, init=0):
        self.data = arr.array(type, [init] * nrows * ncols)
        self.nrows = nrows
        self.ncols = ncols
        self.max = (1 << (self.data.itemsize * 8)) - 1

    def __getitem__(self, pos: "tuple[int,int]"):
        return self.data[pos[0] * self.ncols + pos[1]]

    @property
    def shape(self):
        return self.nrows, self.ncols
    
    @property
    def length(self):
        return len(self.data)
    
    def add_rows(self, nrows: int, init=0) -> int:
        self.data.extend([init] * nrows * self.ncols)
        self.nrows += nrows
        return self.nrows

    def add_set_row(self, values: list):
        if len(values) != self.ncols:
            raise IndexError("length of values doesn't correspond to number of columns in matrix")
        self.data.extend(values)
        self.nrows += 1

    def get_row(self, row: int) -> arr.array:
        if row >= self.nrows:
            raise IndexError("array index out of range")
        start = row * self.ncols
        return self.data[start:(start + self.ncols)]

    def set_row(self, row: int, values):
        if row >= self.nrows:
            raise IndexError("array index out of range")
        start = row * self.ncols
        self.data[start:(start + self.ncols)] = arr.array(self.data.typecode, values)

    def get_col(self, col: int) -> arr.array:
        if col >= self.ncols:
            raise IndexError("array index out of range")
        return arr.array(self.data.typecode, self.data[col::self.ncols])
    
    def set_col(self, col: int, values):
        if col >= self.ncols:
            raise IndexError("array index out of range")
        self.data[col::self.ncols] = arr.array(self.data.typecode, values)

    def get(self, row: int, col: int):
        if row >= self.nrows or col >= self.ncols:
            raise IndexError("array index out of range")
        return self.data[row * self.ncols + col]
    
    def set(self, row: int, col: int, value):
        if row >= self.nrows or col >= self.ncols:
            raise IndexError("array index out of range")
        self.data[row * self.ncols + col] = value


class Depth(Matrix):
    def __init__(self, region: Region, nsamples: int):
        super().__init__('L', region.length, nsamples, 0)
        self.region = region

    def is_valid(self, pos: int, sidx: int, mindepth: int) -> bool:
        return self.get(pos-self.region.start, sidx) >= mindepth

    def get_position(self, pos: int, sidx: int) -> int:
        return self.get(pos-self.region.start, sidx)
    
    def get_position_row(self, pos: int) -> arr.array:
        return self.get_row(pos-self.region.start)

    def set_position(self, pos: int, sidx: int, depth: int):
        value = depth if depth <= self.max else self.max
        self.set(pos-self.region.start, sidx, value)

    def get_colsum(self, spos: int, epos: int, sidx: int) -> int:
        srow = spos-self.region.start
        erow = epos-self.region.start
        start = (srow * self.ncols) + sidx
        end = (erow * self.ncols) + sidx
        return sum(self.data[start:end:self.ncols])


class VariantSite(Matrix):
    def __init__(self, nsamples: int, refbase: Base):
        super().__init__('d', nsamples, 4, 0.)
        self.refbase = refbase
        self.data[self.refbase::4] = arr.array('d', [1.0] * nsamples)
        self.sb = arr.array('L', [0] * nsamples)
        self.sor = arr.array('d', [0.] * nsamples)
        self.dp4 = [(0,0,0,0)] * nsamples

    def update_freqs(self, sidx: int, bases: "tuple[str]", freqs: "tuple[float]"):
        current_freqs = self.get_row(sidx)
        for bstr, freq in zip(bases, freqs):
            base = Base[bstr]
            if base == self.refbase:
                raise ValueError("Bases to update include reference base!")
            if current_freqs[base] > 0.:
                raise ValueError(f"Frequency of alternative allele {base.name} already defined!")
            current_freqs[base] = freq
        current_freqs[self.refbase] = 1. - sum(current_freqs[base] for base in Base if base != self.refbase)
        self.set_row(sidx, current_freqs)

    def get_freqs(self, sidx: int) -> "tuple[Base,float]":
        return ((base, freq) for base, freq in zip(Base, self.get_row(sidx)))
    
    def get_derfreqs(self, sidx: int, ancbase: Base) -> "tuple[Base,float]":
        return ((base, freq) for base, freq in zip(Base, self.get_row(sidx)) \
                             if base != ancbase and freq > 0.)
    
    def get_major_base(self, sidx: int) -> "tuple[Base,float]":
        return max(self.get_freqs(sidx), key=lambda item: item[1])
    
    def is_valid_freq(self, minfreq: float) -> bool:
        max_reffreq = 1. - minfreq
        return any(reffreq <= max_reffreq for reffreq in self.data[self.refbase::4])

    def is_valid_sb(self, maxsb: int) -> bool:
        return any(sb <= maxsb for sb in self.sb)

    def is_valid_sor(self, maxsor: int) -> bool:
        return any(sor <= maxsor for sor in self.sor)
    
    def is_valid_strandcov(self, sidx: int, minscov: int) -> bool:
        return min(self.dp4[sidx][0] + self.dp4[sidx][2], self.dp4[sidx][1] + self.dp4[sidx][3]) >= minscov


class VariantSite_generic(Matrix):
    def __init__(self, chrom: str, pos: int, nsamples: int, refallele: str):
        super().__init__('d', 1, nsamples, 0.)  # start with 1 allele as row for easy extensibility
        self.chrom = chrom
        self.pos = pos
        self.alleles = [refallele]
        self.allele_lookup = {refallele: 0}
        self.data[0:] = arr.array('d', [1.0] * nsamples)
        self.sb = arr.array('L', [0] * nsamples)
        self.sor = arr.array('d', [0.] * nsamples)

    def add_allele(self, allele: str) -> int:
        if allele not in self.allele_lookup:
            aidx = len(self.alleles)   # new allele added, new allele index is equal to length of list before append.
            self.allele_lookup[allele] = aidx
            self.alleles.append(allele)
            self.add_rows(1, 0.)
            return aidx
        else:
            return self.allele_lookup[allele]

    def update_freqs(self, sidx: int, alleles: "tuple[str]", freqs: "tuple[float]"):
        for allele, freq in zip(alleles, freqs):
            aidx = self.add_allele(allele)
            if aidx == 0:
                raise ValueError("Alleles to update include reference allele!")
            if self.get(aidx, sidx) > 0.:
                raise ValueError(f"{self.chrom}:{self.pos} - Frequency of alternative allele {allele} already defined for sample {sidx}!")
            self.set(aidx, sidx, freq)
        ref_freq = 1. - sum(self.get_col(sidx)[1:])
        self.set(0, sidx, ref_freq)

    def get_allele(self, idx: int) -> str:
        if idx >= len(self.alleles):
            raise IndexError("Allele index out of range")
        return self.alleles[idx]
    
    def get_freqs(self, sidx: int):
        return self.get_col(sidx)
    
    def get_derfreqs(self, sidx: int, ancidx: int) -> "tuple[int, float]":
        return ((idx, freq) for idx, freq in enumerate(self.get_col(sidx)) \
                             if idx != ancidx and freq > 0.)
    
    def get_major(self, sidx: int) -> int:
        freqs = self.get_col(sidx)
        return freqs.index(max(freqs))
    
    def get_major_allele(self, sidx: int) -> str:
        return self.alleles[self.get_major(sidx)]
    
    def is_valid_freq(self, minfreq: float) -> bool:
        max_reffreq = 1. - minfreq
        return any(reffreq <= max_reffreq for reffreq in self.get_row(0))

    def is_valid_sb(self, maxsb: int) -> bool:
        return any(sb <= maxsb for sb in self.sb)

    def is_valid_sor(self, maxsor: int) -> bool:
        return any(sor <= maxsor for sor in self.sor)


class Diversity(Matrix):
    def __init__(self, region: Region, nsamples: int):
        super().__init__('d', nsamples, 6, 0.)
        self.nvalid = arr.array('H', [0] * nsamples)
        self.region = region
        self.stats = {'pi1': 0, 'pi2': 1, 'pitot': 2, 'dxy': 3, 'h':4, 'fst': 5}

    def calculate_pi(self, freqs1: Tuple[float, ...], freqs2: Tuple[float, ...]) -> Tuple[float, ...]:
        pi1 = 1 - sum(freq * freq for freq in freqs1)
        pi2 = 1 - sum(freq * freq for freq in freqs2)
        pitot = 1 - sum(((freq1 + freq2)/2)**2 for freq1, freq2 in zip(freqs1, freqs2))
        dxy = 1 - sum(freq1 * freq2 for freq1, freq2 in zip(freqs1, freqs2))
        h = -sum(freq * math.log(freq) for freq in freqs2 if freq > 0.)
        return pi1, pi2, pitot, dxy, h

    def add_site(self, sidx: int, ancfreqs: Tuple[float, ...], freqs: Tuple[float, ...]):
        pis = self.calculate_pi(ancfreqs, freqs)
        for vidx, val in enumerate(pis):
            self.data[sidx * self.ncols + vidx] += val
        self.nvalid[sidx] += 1

    def get_diversity(self, minprop: float=0.1) -> Matrix:
        minvalid = self.region.length * minprop
        divstats = Matrix('d', self.nrows, self.ncols, 0.)
        for sidx in range(self.nrows):
            if self.nvalid[sidx] >= minvalid:
                stats = self.get_row(sidx)
                corr = [stat / self.nvalid[sidx] for stat in stats[:-1]]
                fst = (corr[2] - (corr[0] + corr[1])/2) / corr[2] if corr[2] else float('nan')
                corr.append(fst)
                divstats.set_row(sidx, corr)
            else:
                divstats.set_row(sidx, [float('nan')] * self.ncols)
        return divstats


@dataclass
class Variant:
    ref: str
    major: str
    depth: int
    freqs: "tuple[float,float,float,float]"
    bases: "tuple[tuple]"
    sb: int
    dp4: "tuple[int,int,int,int]"


def get_sor(ref_fwd: int, ref_rev: int, alt_fwd:int, alt_rev: int) -> float:
    ref_fwd += 1
    ref_rev += 1
    alt_fwd += 1
    alt_rev += 1
    R = (ref_fwd * alt_rev) / (alt_fwd * ref_rev)
    sym_ratio = R + (1 / R)
    ref_ratio = min(ref_fwd, ref_rev) / max(ref_fwd, ref_rev)
    alt_ratio = min(alt_fwd, alt_rev) / max(alt_fwd, alt_rev)
    return math.log(sym_ratio) + math.log(ref_ratio) - math.log(alt_ratio)


def get_sampleids(vcf_files: "list[Path]") -> "list[str]":
    return list(fn.name.split('.')[0] for fn in vcf_files)


def get_sample_info(samples: "list[str]") -> "dict[str,tuple[str,int,str]]":
    sampleinfo = {}
    pattern = re.compile(r'^(\w+)_P(\d+)_([A-Z0-9](_new)?)$')
    for sample in samples:
        m = pattern.match(sample)
        try:
            sampleid = m.group(1)
            passage = int(m.group(2))
            replicate = m.group(3)
        except AttributeError:
            logger.warning(f"Sample string {sample} doesn't have a valid structure.")
            fields = sample.split('_')
            sampleid = fields[0]
            passage = -1
            replicate = 'A'
        sampleinfo[sample] = (sampleid, passage, replicate)
    return sampleinfo


def read_fasta(fasta_file: Path) -> dict:
    with FastaFile(str(fasta_file)) as fasta_in:
        seqs = {ref: fasta_in.fetch(ref) for ref in fasta_in.references}
    return seqs

def read_fai(fai_file: Path) -> "list[Region]":
    regions = []
    with open(fai_file, 'r') as fai_in:
        for line in fai_in:
            fields = line.strip().split('\t')
            if not len(fields) == 5:
                raise Exception("Invalid format of line in fasta index!")
            regions.append(Region(fields[0], 1, int(fields[1])))
    return regions


def get_sequence(fasta_file: Path, region: Region) -> str:
    with FastaFile(str(fasta_file)) as fasta_in:
        seq = fasta_in.fetch(region.chrom, region.start-1, region.end)
    return seq


def get_contigs(tabix_file: Path) -> list:
    with VariantFile(str(tabix_file)) as vcf_in:
        return vcf_in.contigs


def get_samples_from_depth_file(depth_file: Path) -> "list[str]":
    with TabixFile(str(depth_file)) as depth_in:
        header = depth_in.header[0].split('\t')
    return [Path(fn).name.split('.')[0] for fn in header[2:]]


def get_column_dict(header: "list[str]", samples: "list[str]", validate: bool=True) -> dict:
    columns_by_id = {Path(fn).name.split('.')[0]: col for col, fn in enumerate(header) if col >= 2}
    for sample in samples:
        if sample in columns_by_id:
            logger.info(f"Found {sample} at column {columns_by_id[sample] + 1} in depth file.")
        elif validate:
            raise Exception(f"Sample {sample} not found in depth file!")
        else:
            logger.warning(f"Sample {sample} not found in depth file!")
    return columns_by_id


def get_ancestral_sequence(
        refseq: str,
        region: Region,
        depths: Depth,
        variants: "dict[int,VariantSite]",
        anc_idx: int,
        mindepth: int=100
        ) -> str:
    if not region.defined:
        raise ValueError("get_ancestral_sequence requires a valid region!")
    if not len(refseq) == region.length:
        raise ValueError(f"reference sequence length {len(refseq)} doesn't match provided region {region.length}!")
    ancseq = arr.array('u', 'N' * region.length)
    for offset, refbase in enumerate(refseq):
        pos = region.start + offset
        if depths.get(offset, anc_idx) < mindepth: continue
        if pos in variants:
            if refbase.upper() != variants[pos].refbase.name:
                logger.warning("Reference base in reference sequence doesn't match reference base in vcf file.")
            ancseq[offset] = variants[pos].get_major_base(anc_idx)[0].name
        else:
            ancseq[offset] = refbase
    return ancseq.tounicode()


def get_ancestral_sequence_generic(
        refseq: str,
        region: Region,
        depths: Depth,
        variants: "dict[int,VariantSite_generic]",
        anc_idx: int=None,
        mindepth: int=100
        ) -> str:
    if not region.defined:
        raise ValueError("get_ancestral_sequence requires a valid region!")
    if not len(refseq) == region.length:
        raise ValueError(f"reference sequence length {len(refseq)} doesn't match provided region {region.length}!")
    ancseq = [''] * region.length  # This is inefficient, but elements need to be able to contain arbitrary long insertion alleles.
    refmod = 0
    for offset, refbase in enumerate(refseq):
        pos = region.start + offset
        refmod -= 1 # if the reference allele is longer than 1, it will effectively overwrite also the next bases.
        valid = False if anc_idx is not None and depths.get(offset, anc_idx) < mindepth else True
        if pos in variants:
            if refmod > 0:
                logger.warning(f"Overlapping variants detected at position {pos}!")
                continue
            ref_allele = variants[pos].alleles[0]
            major_allele = variants[pos].get_major_allele(anc_idx) if anc_idx is not None else ref_allele
            max_length = max(len(allele) for allele in variants[pos].alleles)
            if len(ref_allele) > 1:
                refmod = len(ref_allele)
            final_allele = major_allele + '-' * (max_length - len(major_allele)) if valid else 'N' * max_length
            ancseq[offset] = final_allele
        else:
            if refmod > 0:  # previous reference allele is still valid for this position
                continue
            else:
                ancseq[offset] = refbase if valid else 'N'
    return "".join(ancseq)


def read_depth_array_list(
        depth_files: "list[Path]",
        region: Region,
        samples: "list[str]",
        interval: int=1000
        ) -> Matrix:
    if not region.defined:
        raise ValueError("read_depth_array requires a valid region!")
    depths = Depth(region, len(samples))
    for depth_file in depth_files:
        with TabixFile(str(depth_file)) as depth_in:
            try:
                header = depth_in.header[0].split('\t')
                columns_by_id = get_column_dict(header, samples, validate=False)
            except Exception as e:
                logger.error(e)
                sys.exit(1)
            processed = 0
            for row in depth_in.fetch(region.chrom, region.start-1, region.end, parser=asTuple()):
                pos = int(row[1])
                for sidx, sample in enumerate(samples):
                    if sample in columns_by_id:
                        depth = int(row[columns_by_id[sample]])
                        depths.set_position(pos, sidx, depth)
                processed += 1
                if not processed % interval: logger.info(f"Processed {processed} lines of depth file.")
            logger.info(f"Processed {processed} lines of depth file.")
    return depths


def read_depth_array(
        depth_file: Path,
        region: Region,
        samples: "list[str]",
        interval: int=1000
        ) -> Matrix:
    if not region.defined:
        raise ValueError("read_depth_array requires a valid region!")
    depths = Depth(region, len(samples))
    with TabixFile(str(depth_file)) as depth_in:
        try:
            header = depth_in.header[0].split('\t')
            columns_by_id = get_column_dict(header, samples)
        except Exception as e:
            logger.error(e)
            sys.exit(1)
        processed = 0
        for row in depth_in.fetch(region.chrom, region.start-1, region.end, parser=asTuple()):
            pos = int(row[1])
            for sidx, sample in enumerate(samples):
                depth = int(row[columns_by_id[sample]])
                depths.set_position(pos, sidx, depth)
            processed += 1
            if not processed % interval: logger.info(f"Processed {processed} lines of depth file.")
        logger.info(f"Processed {processed} lines of depth file.")
    return depths


def read_depth(depth_file: Path, region: Region, sample_id: str, seqs=None, interval: int=1000):
    if not seqs is None:
        depths = {seqname: [0] * len(seq) for seqname, seq in seqs.items()}
    else:
        depths = {}
    with TabixFile(str(depth_file)) as depth_in:
        try:
            header = depth_in.header[0].split('\t')
            columns_by_id = get_column_dict(header, [sample_id])
        except Exception as e:
            logger.error(e)
            sys.exit(1)
        samplecol = columns_by_id[sample_id]
        processed = 0
        for row in depth_in.fetch(region.chrom, region.start-1, region.end, parser=asTuple()):
            pos = int(row[1])
            if row[0] not in depths:
                depths[row[0]] = [0] * 1000
            while len(depths[row[0]]) < pos:
                depths[row[0]].extend([0] * 1000)
            depths[row[0]][pos-1] = int(row[samplecol])
            processed += 1
            if not processed % interval: logger.info(f"Processed {processed} lines.")
        logger.info(f"Processed {processed} lines of depth file.")
    return depths


# Under construction!
def read_multisample_vcf(
        vcf_file: Path,
        region: Region,
        samples: "list[str]",
        interval: int=1000
        ) -> "dict[VariantSite_generic]":
    if not region.defined:
        raise ValueError("read_multisample_vcf requires a valid region!")
    variants = {}
    with VariantFile(str(vcf_file)) as vcf_in:
        processed = 0
        for rec in vcf_in.fetch(region.chrom, region.start-1, region.end):
            pos = int(rec.pos)
            variants[pos] = VariantSite_generic(region.chrom, pos, len(samples), rec.ref)   # Need a different data type here.
            for sidx, sample in enumerate(samples):
                gt = rec.samples[sample]["GT"]
                dp = rec.samples[sample]["DP"]
                variants[pos].add_genotype(sidx, rec.alts, gt)
            processed += 1
            if not processed % interval: logger.info(f"Processed {processed} lines.")
        logger.info(f"Processed {processed} lines of VCF file {vcf_file}.")
    return variants


# LoFreq doesn't adhere to VCF standard and puts multiallelic variants on multiple lines in the vcf file! 
def read_vcf_array(
        vcf_files: "list[Path]",
        samples: "list[str]",
        region: Region,
        max_sb: int=100,
        interval: int=1000
        ) -> "dict[VariantSite]":
    if not region.defined:
        raise ValueError("read_vcf_array requires a valid region!")
    variants = {}
    for sidx, vcf_file in enumerate(vcf_files):
        with VariantFile(str(vcf_file)) as vcf_in:
            # Since LoFreq doesn't produce proper VCF headers, fetch with contig will fail if the vcf file is empty.
            if vcf_in.get_tid(region.chrom) < 0: continue
            processed = 0
            for rec in vcf_in.fetch(region.chrom, region.start-1, region.end):
#               Filtering for strand bias like this is a bad idea, as it just sets the frequency of any dubious variant to zero.
#                if int(rec.info["SB"]) > max_sb: continue
                pos = int(rec.pos)
                refbase = Base[rec.ref]
                if pos not in variants:
                    variants[pos] = VariantSite(len(vcf_files), refbase)
                elif refbase != variants[pos].refbase:
                    logger.critical(f"Inconsistent reference alleles ({refbase.name} vs. {variants[pos].refbase.name}) \
                                    for variant {rec.chrom}:{pos}!")
                altfreqs = tuple(float(freq) for freq in str(rec.info["AF"]).split(','))
                variants[pos].update_freqs(sidx, rec.alts, altfreqs)
                if int(rec.info["SB"]) > variants[pos].sb[sidx]:
                    variants[pos].sb[sidx] = int(rec.info["SB"])
                variants[pos].dp4[sidx] = rec.info["DP4"]
                sor = get_sor(*rec.info["DP4"])
                if sor > variants[pos].sor[sidx]:
                    variants[pos].sor[sidx] = sor
                processed += 1
                if not processed % interval: logger.info(f"Processed {processed} lines.")
            logger.info(f"Processed {processed} lines of VCF file {vcf_file}.")
    return variants


# LoFreq doesn't adhere to VCF standard and puts multiallelic variants on multiple lines in the vcf file! 
def read_vcf_array_generic(
        vcf_files: "list[Path]",
        region: Region,
        interval: int=1000
        ) -> "dict[VariantSite_generic]":
    if not region.defined:
        raise ValueError("read_vcf_array_generic requires a valid region!")
    variants = {}
    for sidx, vcf_file in enumerate(vcf_files):
        with VariantFile(str(vcf_file)) as vcf_in:
            # Since LoFreq doesn't produce proper VCF headers, fetch with contig will fail if the vcf file is empty.
            if vcf_in.get_tid(region.chrom) < 0: continue
            processed = 0
            for rec in vcf_in.fetch(region.chrom, region.start-1, region.end):
                pos = int(rec.pos)
                if pos not in variants:
                    variants[pos] = VariantSite_generic(region.chrom, pos, len(vcf_files), rec.ref)
                elif rec.ref != variants[pos].alleles[0]:
                    logger.critical(f"Inconsistent reference alleles ({rec.ref} vs. {variants[pos].alleles[0]}) for variant {rec.chrom}:{pos}!")
                    continue
                altfreqs = tuple(float(freq) for freq in str(rec.info["AF"]).split(','))
                variants[pos].update_freqs(sidx, rec.alts, altfreqs)
                if int(rec.info["SB"]) > variants[pos].sb[sidx]:
                    variants[pos].sb[sidx] = int(rec.info["SB"])
                processed += 1
                if not processed % interval: logger.info(f"Processed {processed} lines.")
            logger.info(f"Processed {processed} lines of VCF file {vcf_file}.")
    return variants


# LoFreq doesn't adhere to VCF standard and puts multiallelic variants on multiple lines in the vcf file! 
def read_vcf(vcf_file: Path, region: Region, sample_id: str, interval: int=1000):
    variants = {}
    with VariantFile(str(vcf_file)) as vcf_in:
        # Since LoFreq doesn't produce proper VCF headers, fetch with contig will fail if the vcf file is empty.
        if region.chrom is not None and vcf_in.get_tid(region.chrom) < 0:
            return variants
        processed = 0
        for rec in vcf_in.fetch(region.chrom, region.start-1, region.end):
            if rec.chrom not in variants:
                variants[rec.chrom] = {}
            depth = int(rec.info["DP"])
            if rec.pos in variants[rec.chrom]:
                variant = variants[rec.chrom][rec.pos]
                freqs_by_base = dict(zip(BASES, variant.freqs))
                if rec.ref != variant.ref:
                    logger.warning(f"Inconsistent reference alleles ({rec.ref} vs {variant.ref}) \
                                     for {rec.chrom}:{rec.pos} despite definition on multiple lines!")
                if len(rec.alts) > 1:
                    logger.warning(f"Multiple alternative allele entries \
                                     for {rec.chrom}:{rec.pos} despite definition on multiple lines!")
                altbase = rec.alts[0]
                if freqs_by_base[altbase]:
                    logger.warning(f"Frequency of alternative allele {altbase} already defined \
                                     for {rec.chrom}:{rec.pos} despite definition on multiple lines!")
                freqs_by_base[altbase] = rec.info["AF"]
                freqs_by_base[rec.ref] = 1 - sum(freq for base, freq in freqs_by_base.items() if base != rec.ref)
            else:
                if isinstance(rec.info["AF"], float):
                    freqs = (1 - rec.info["AF"], rec.info["AF"])
                else:
                    alt_freqs = tuple(float(freq) for freq in rec.info["AF"].split(','))
                    freqs = (1 - sum(alt_freqs), *alt_freqs)
                freqs_by_base = {base: freq for base, freq in zip(rec.alleles, freqs)}
            freqs_ordered = tuple(freqs_by_base.get(base, 0.) for base in BASES)
            sorted_bases = [(base, freq) for base, freq in sorted(freqs_by_base.items(), key=lambda item: item[1], reverse=True)]
            variants[rec.chrom][rec.pos] = Variant(rec.ref, sorted_bases[0][0], depth, freqs_ordered,
                                                   sorted_bases, int(rec.info["SB"]), rec.info["DP4"])
#            variants[rec.chrom][rec.pos] = (rec.ref, sorted_bases[0][0], depth, freqs_ordered, sorted_bases, int(rec.info["SB"]), rec.info["DP4"])
            processed += 1
            if not processed % interval: logger.info(f"Processed {processed} lines.")
        logger.info(f"Processed {processed} lines of VCF file.")
    return variants


def read_gtf(gtf_file: Path, interval: int=1000) -> "list[Region]":
    regions = []
    with open(gtf_file, 'r') as inhandle:
        idpattern = re.compile(r'gene_id "(.+)"')
        processed = 0
        for line in inhandle:
            fields = line.strip().split('\t')
            if len(fields) != 9:
                logger.error(f"Wrongly formatted line in GTF file!")
                sys.exit(1)
            m = idpattern.match(fields[8])
            try:
                region_name = m.group(1)
            except AttributeError:
                logger.error(f"Feature doesn't have a valid tag ({fields[8]}). Skipping line.")
                continue
            regions.append(Region(fields[0], int(fields[3]), int(fields[4]), region_name))
            processed += 1
            if not processed % interval: logger.info(f"Processed {processed} lines.")
        logger.info(f"Processed {processed} lines of GTF file.")
    return regions


def read_bed(bed_file: Path, interval: int=1000) -> "list[Region]":
    regions = []
    with open(bed_file, 'r') as inhandle:
        processed = 0
        for line in inhandle:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                logger.error(f"Wrongly formatted line in BED file!")
                sys.exit(1)
            regions.append(Region(fields[0], int(fields[1]) + 1, int(fields[2])))
            processed += 1
            if not processed % interval: logger.info(f"Processed {processed} lines.")
        logger.info(f"Processed {processed} lines of BED file.")
    return regions


def read_samplelist(sample_file: Path) -> "dict[str, tuple]":
    samples = {}
    with open(sample_file, 'r') as inhandle:
        for line in inhandle:
            fields = line.strip().split('\t')
            if len(fields) != 5:
                logger.error(f"Wrongly formatted line in sample list file!")
                sys.exit(1)
            samples[fields[0]] = {'setting': fields[1], 'passage': int(fields[2]), 'replicate': fields[3], 'vcf': Path(fields[4])}
    return samples


def get_regions(bed_file: Path=None, gtf_file: Path=None, fai_file: Path=None, region_strings: "list[str]"=None) -> "list[Region]":
    regions = []
    if region_strings:
        for region_str in region_strings:
            try:
                region = Region.from_string(region_str)
            except ValueError as e:
                logger.critical(e)
                sys.exit(1)
            regions.append(region)
    if bed_file:
        regions.extend(read_bed(bed_file))
    if gtf_file:
        regions.extend(read_gtf(gtf_file))
    if fai_file:
        regions.extend(read_fai(fai_file))
    if len(regions) == 0:
        logger.error(f"No valid region provided!")
        sys.exit(1)
    return regions
