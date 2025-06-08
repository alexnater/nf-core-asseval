#!/usr/bin/env python

import sys
import argparse
import logging
import csv
import array as arr
from pathlib import Path
from collections import Counter
from pysam import VariantFile
from divlib import Region, Depth
from divlib import get_sample_info, get_sequence, get_samples_from_depth_file, read_depth_array, get_regions

logger = logging.getLogger()


class VariantSite:
    def __init__(
            self,
            chrom: str,
            pos: int,
            qual: float,
            samples: list[str],
            depths: list[int],
            gqs: list[int],
            alleles: list[int],
            phased: list[int]):
        self.chrom = chrom
        self.pos = pos
        self.qual = qual
        self.depths = arr.array('L', depths) if depths is not None else arr.array('H', [0] * len(samples))
        self.gqs = arr.array('H', gqs) if gqs is not None else arr.array('H', [0] * len(samples))
        self.alleles = arr.array('b', alleles) if alleles is not None else arr.array('b', [-1] * 2 * len(samples))
        self.phased = arr.array('b', phased) if phased is not None else arr.array('b', [False] * len(samples))

    def set_individual(
            self,
            sidx: int,
            qual: int,
            depth: int,
            gq: int,
            alleles: tuple[int,int],
            phased: bool):
        if sidx >= len(self.depths):
            raise Exception("Sample index out of bounds")
        self.qual = qual if self.qual is None or qual < self.qual else self.qual
        self.depths[sidx] = depth
        self.gqs[sidx] = gq
        self.alleles[2*sidx] = alleles[0]
        self.alleles[2*sidx+1] = alleles[1]
        self.phased[sidx] = phased


def calculate_window_stats(
        region: Region,
        depths: Depth,
        variants: dict[str,VariantSite],
        refseq: str,
        windowsize: int,
        stepsize: int,
        mindepth: int,
        mingq: int,
        total_dp: list[float]
        ) -> "tuple[arr.array,arr.array,arr.array,arr.array]":
    cgcs, cats = 0., 0.
    winvalues = []
    winstart = region.start - 1
    while winstart < region.end:
        winend = winstart + windowsize
        winend = winend if winend < region.end else region.end
        if not refseq is None:
            seqstart = winstart - region.start + 1
            seqend = winend - region.start + 1
            counter = Counter(refseq[seqstart:seqend].upper())
            totvalid = (counter['A'] + counter['T'] + counter['G'] + counter['C'])
            gc = (counter['G'] + counter['C']) / totvalid if totvalid > 0 else float('nan')
            if (counter['G'] + counter['C']) > 0:
                gc_skew = (counter['G'] - counter['C']) / (counter['G'] + counter['C'])
                cgcs += gc_skew
            if (counter['A'] + counter['T']) > 0:
                at_skew = (counter['A'] - counter['T']) / (counter['A'] + counter['T'])
                cats += at_skew
        else:
            gc = float('nan')
            gc_skew, at_skew = float('nan'), float('nan')
        nvalid = arr.array('L', [0] * depths.ncols)
        covsums = arr.array('L', [0] * depths.ncols)
        hetsums = arr.array('L', [0] * depths.ncols)
        logger.info(f"Working on window {region.chrom}:{winstart+1}-{winend} ...")
        for pos in range(winstart + 1, winend + 1):
            valid_samples = [depth >= mindepth for depth in depths.get_position_row(pos)]
            for sidx, is_valid in enumerate(valid_samples):
                if is_valid:
                    nvalid[sidx] += 1
            if not pos in variants:
                continue
            alleles = (al for al in variants[pos].alleles)
            hets = [ True if (is_valid and gq >= mingq and al1 >= 0 and al2 >= 0 and al1 != al2) else False \
                    for al1, al2, gq, is_valid in zip(alleles, alleles, variants[pos].gqs, valid_samples) ]
            for sidx, is_het in enumerate(hets):
                if is_het:
                    hetsums[sidx] += 1
        for sidx in range(depths.ncols):
            covsums[sidx] = depths.get_colsum(winstart + 1, winend + 1, sidx)
        winvalues.append([region.chrom, winstart, winend, gc, gc_skew, cgcs, at_skew, cats, \
                          nvalid, hetsums, \
                          (covsum / windowsize for covsum in covsums), \
                          (covsum / (windowsize * total_dp[sidx]) for sidx, covsum in enumerate(covsums))])
        winstart += stepsize
    return winvalues


def read_vcf(
        vcf_file: Path,
        region: Region,
        samples: list[str],
        include_indels: bool=True,
        interval: int=1000
        ):
    variants = {}
    with VariantFile(str(vcf_file)) as vcf_in:
        samples = list(vcf_in.header.samples) if samples is None else samples
        if vcf_in.get_tid(region.chrom) < 0:
            return samples, variants
        processed = 0
        for rec in vcf_in.fetch(region.chrom, region.start-1, region.end):
            depths = [rec.samples[sample]["DP"] if rec.samples[sample]["DP"] is not None else 0 for sample in samples]
            gqs = [rec.samples[sample]["GQ"] if rec.samples[sample]["GQ"] is not None else 0 for sample in samples]
            alleles = [al if al is not None else -1 for sample in samples for al in rec.samples[sample]["GT"]]
            phased = [rec.samples[sample].phased for sample in samples]
            if rec.rlen == 1 or (rec.rlen > 1 and include_indels):
                variants[rec.pos] = VariantSite(rec.chrom, rec.pos, rec.qual, samples, depths, gqs, alleles, phased)
            processed += 1
            if not processed % interval: logger.info(f"Processed {processed} lines.")
        logger.info(f"Processed {processed} lines of VCF file {vcf_file}.")
    return samples, variants

def read_vcf_single(
        vcf_files: list[Path],
        region: Region,
        samples: list[str],
        include_indels: bool=True,
        interval: int=1000
        ):
    sample_dict = {sample: sidx for sidx, sample in enumerate(samples)}
    variants = {}
    for vcf_file in vcf_files:
        with VariantFile(str(vcf_file)) as vcf_in:
            sample = list(vcf_in.header.samples)[0]
            try:
                sidx = sample_dict[sample]
            except:
                raise Exception(f"VCF sample id {sample} is not in list of sample names")
            if vcf_in.get_tid(region.chrom) < 0: continue
            processed = 0
            for rec in vcf_in.fetch(region.chrom, region.start-1, region.end):
                depth = rec.samples[sample]["DP"] if rec.samples[sample]["DP"] is not None else 0
                gq = rec.samples[sample]["GQ"] if rec.samples[sample]["GQ"] is not None else 0
                alleles = tuple(al if al is not None else -1 for al in rec.samples[sample]["GT"])
                phased = rec.samples[sample].phased
                if rec.rlen == 1 or (rec.rlen > 1 and include_indels):
                    if not rec.pos in variants:
                        variants[rec.pos] = VariantSite(rec.chrom, rec.pos, None, samples, None, None, None, None)
                    variants[rec.pos].set_individual(sidx, rec.qual, depth, gq, alleles, phased)
                processed += 1
                if not processed % interval: logger.info(f"Processed {processed} lines.")
            logger.info(f"Processed {processed} lines of VCF file {vcf_file}.")
    return samples, variants


def read_mosdepth(infiles: list[Path]) -> tuple[dict, list, dict]:
    total = [0.] * len(infiles)
    depth = {}
    chrlen = {}
    for fidx, infile in enumerate(infiles):
        with open(infile, 'r', newline='') as inhandle:
            reader = csv.DictReader(inhandle, delimiter='\t')
            for row in reader:
                if row['chrom'] == 'total':
                    total[fidx] = float(row['mean'])
                if not row['chrom'] in chrlen:
                    chrlen[row['chrom']] = int(row['length'])
                elif row['chrom'] != 'total' and int(row['length']) != chrlen[row['chrom']]:
                    logger.error(f"ERROR: Unequal length of reference contig {row['chrom']}: {row['length']} vs. {chrlen[row['chrom']]}")
                if not row['chrom'] in depth:
                    depth[row['chrom']] = [0.] * len(infiles)
                depth[row['chrom']][fidx] = float(row['mean'])
    return depth, total, chrlen


def print_summary(outfile: Path, depth: dict, total: list[float], chrom_lengths: dict, samples: list[str]):
    with open(outfile, 'w') as outhandle:
        print("contig", "length", "\t".join(f"raw_{sample.split('_')[0]}" for sample in samples), \
              "\t".join(f"corr_{sample.split('_')[0]}" for sample in samples), "F/M", sep="\t", file=outhandle)
        for contig, dps in depth.items():
            print(contig, chrom_lengths[contig], "\t".join(f"{dp:.2f}" for dp in dps), \
                  "\t".join(f"{dp / tot:.2f}" for dp, tot in zip(dps, total)), \
                    f"{(dps[2] / dps[1]) * (total[1] / total[2]):.2f}" if dps[1] > 0 else float("nan"), sep="\t", file=outhandle)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculates windowwise corrected mean depth per sample by region.",
        epilog="Example: python windows_stats.py depth.tsv out1 --samples sample1 sample2 sample3 --bed regions.bed --windowsize 10000 --stepsize 2000",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--depth",
        type=Path,
        required=True,
        help="Site-wise coverage information in bgzipped and tabixed SAMtools depth format."
    )
    parser.add_argument(
        "--vcf",
        type=Path,
        nargs='+',
        required=True,
        help="VCF file with joint variant calls or list of single-sample VCF files, bgzipped and tabixed."
    )
    parser.add_argument(
        "-o",
        "--outprefix",
        type=str,
        required=True,
        help="Prefix for output files."
    )
    parser.add_argument(
        "-s",
        "--samples",
        nargs='+',
        type=str,
        required=False,
        help="List of sample IDs."
    )
    parser.add_argument(
        "-r",
        "--regions",
        nargs='+',
        type=str,
        required=False,
        help="List of regions in chrom:start-end format to work on.",
    )
    parser.add_argument(
        "--summaries",
        nargs='+',
        type=Path,
        required=False,
        help="List of mosdepth summary files.",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=False,
        help="Reference genome FASTA file."
    )
    parser.add_argument(
        "--fai",
        type=Path,
        required=False,
        help="Fasta index of reference genome."
    )
    parser.add_argument(
        "--bed",
        type=Path,
        required=False,
        help="BED file with regions to process."
    )
    parser.add_argument(
        "--gtf",
        type=Path,
        required=False,
        help="GTF file with regions to process."
    )
    parser.add_argument(
        "--minlength",
        type=int,
        help="Minium length of region to consider.",
        default=1000000
    )
    parser.add_argument(
        "-d",
        "--mindepth",
        type=int,
        help="Minium depth needed to be considered as valid genotype call.",
        default=5
    )
    parser.add_argument(
        "-q",
        "--mingq",
        type=int,
        help="Minium genotype quality needed to be considered as valid genotype call.",
        default=0
    )
    parser.add_argument(
        "--include_indels",
        help="Include INDEL variants.",
        action='store_true'
    )
    parser.add_argument(
        "-w",
        "--windowsize",
        type=int,
        help="Size of windows.",
        default=100
    )
    parser.add_argument(
        "--stepsize",
        type=int,
        help="Stepsize for sliding windows.",
        default=100
    )
    parser.add_argument(
        "--total_depth",
        nargs='+',
        type=float,
        required=False,
        help="Total genomic mean depth for each sample.",
    )
    parser.add_argument(
        "-i",
        "--interval",
        type=int,
        help="The desired interval of progress updates.",
        default=1000
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING"
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    regions = [region for region in get_regions(args.bed, args.gtf, args.fai, args.regions) if region.length >= args.minlength]
    samples = args.samples if args.samples else get_samples_from_depth_file(args.depth)
    total_depth = args.total_depth if args.total_depth else [1.] * len(samples)

    if args.summaries:
        chrom_depth, mean_dp, chrom_lengths = read_mosdepth(args.summaries)
        print_summary(Path(f"{args.outprefix}.summary.tsv"), chrom_depth, mean_dp, chrom_lengths, samples)
        if not args.total_depth:
            logger.info(f"Setting mean depth to {mean_dp}.")
            total_depth = mean_dp

    with open(Path(f"{args.outprefix}_win{args.windowsize}_step{args.stepsize}.bed"), 'w') as outhandle:
        print("#chrom", "start", "end", "gc_content", "gc_skew", "cgc_skew", "at_skew", "cat_skew", \
            "\t".join(f"{sample.split('_')[0]}_val" for sample in samples), \
            "\t".join(f"{sample.split('_')[0]}_het" for sample in samples), \
            "\t".join(f"{sample.split('_')[0]}_raw" for sample in samples), \
            "\t".join(f"{sample.split('_')[0]}_corr" for sample in samples), sep="\t", file=outhandle)

        for region in regions:
            if args.fasta:
                try:
                    refseq = get_sequence(args.fasta, region)
                except KeyError as e:
                    logger.critical(f"Couldn't find {region} in reference sequence!")
                    sys.exit(1)
            else:
                refseq = None
            
            logger.info(f"Working on region {region} ...")
            depths = read_depth_array(args.depth, region, samples, args.interval)
            if len(args.vcf) == 1:
                _, variants = read_vcf(args.vcf, region, samples, args.include_indels, args.interval)
            else:
                _, variants = read_vcf_single(args.vcf, region, samples, args.include_indels, args.interval)
            winvalues = calculate_window_stats(region, depths, variants, refseq, args.windowsize, \
                                            args.stepsize, args.mindepth, args.mingq, total_depth)
            for chrom, start, end, gc, gc_skew, cgcs, at_skew, cats, valid, hets, dps_raw, dps_corr in winvalues:
                print(chrom, str(start), str(end), f"{gc:.2f}", f"{gc_skew:.2f}", f"{cgcs:.2f}", f"{at_skew:.2f}", f"{cats:.2f}", \
                    "\t".join(f"{val}" for val in valid), "\t".join(f"{het}" for het in hets), \
                    "\t".join(f"{cov:.4f}" for cov in dps_raw), \
                    "\t".join(f"{cov:.4f}" for cov in dps_corr), sep="\t", file=outhandle)


if __name__ == "__main__":
    sys.exit(main())
