#!/usr/bin/env python

import sys
import argparse
from pathlib import Path
import logging
import yaml
import csv
import re

logger = logging.getLogger()


def get_infiles(flagstats: "list[Path]", summaries: "list[Path]", coverages: "list[Path]") -> "dict[str, list[Path]]":
    infiles = {}
    samples = []
    for fn_list in flagstats, summaries, coverages:
        for fn in fn_list:
            sample = fn.name.split('.')[0]
            samples.append(sample)
            if sample in infiles:
                infiles[sample].append(fn)
            else:
                infiles[sample] = [fn]
    return infiles, samples


def check_infiles(infiles: dict):
    for sample, fn_list in infiles.items():
        if not len(fn_list) == 3:
            raise Exception(f"Incomplete list of infiles for sample {sample}!")


def get_counts(samples: "list[str]", before: "list[int]", after: "list[int]") -> "dict[str, tuple[int]]":
    if not len(samples) == len(before) and len(samples) == len(after):
        raise Exception(f"Length of counts doesn't match number of samples!")
    return { sample: (bc, ac) for sample, bc, ac in zip(samples, before, after) }


def read_flagstats(infile: Path) -> dict:
    flagstats = {}
    with open(infile, 'r') as inhandle:
        for line in inhandle:
            m = re.match('^(\d+).+primary$', line)
            if m:
                filtered = int(m.group(1))
                continue
            m = re.match('^(\d+).+primary duplicates$', line)
            if m:
                dups = int(m.group(1))
                continue
            m = re.match('^(\d+).+primary mapped', line)
            if m:
                mapped = int(m.group(1))
                continue
            m = re.match('^(\d+).+properly paired', line)
            if m:
                proppaired = int(m.group(1))
                continue
        
    flagstats['%Filtered'] = "N/A" if filtered == "N/A" or filtered == 0 else f"{(mapped / filtered):.4f}"
    flagstats['%Mapped'] = "N/A" if filtered == "N/A" or filtered == 0 else f"{(mapped / filtered):.4f}"
    flagstats['%ProperlyPaired'] = "N/A" if mapped == "N/A" or mapped == 0 else f"{(proppaired / mapped):.4f}"
    flagstats['%Duplicates'] = "N/A" if mapped == "N/A" or mapped == 0 else f"{(dups / mapped):.4f}"
    flagstats['FilteredReads'] = filtered
    flagstats['Duplicates'] = dups
    flagstats['MappedReads'] = mapped
    flagstats['ProperlyPaired'] = proppaired
    return flagstats


def read_mosdepth(infile: Path) -> dict:
    depth = {}
    with open(infile, 'r', newline='') as inhandle:
        reader = csv.DictReader(inhandle, delimiter='\t')
        for row in reader:
            if row['chrom'] == 'total_region':
                depth['MeanDepth'] = f"{float(row['mean']):.2f}"
                depth['MinDepth'] = int(row['min'])
                depth['MaxDepth'] = int(row['max'])
                break
    return depth


def read_coverage(infile: Path) -> dict:
    coverages = [0] * len(DEPTH_TRHESHOLDS)
    with open(infile, 'r', newline='') as inhandle:
        reader = csv.DictReader(inhandle, fieldnames=('chrom', 'depth', 'coverage'), delimiter='\t')
        prevcov = 0.
        i = 0
        for row in reader:
            if row['chrom'] == 'total':
                depth = int(row['depth'])
                cov = float(row['coverage'])
                if depth < DEPTH_TRHESHOLDS[i]:
                    coverages[i] = prevcov
                    i += 1
                elif depth == DEPTH_TRHESHOLDS[i]:
                    coverages[i] = cov
                    i += 1
                prevcov = cov
                if i == len(DEPTH_TRHESHOLDS):
                    break
            if i < len(DEPTH_TRHESHOLDS):
                coverages[i] = prevcov
    covstats = {f"Coverage{DEPTH_TRHESHOLDS[i]}": f"{coverage:.2f}" for i, coverage in enumerate(coverages)}
    return covstats


def read_stats(infiles: dict, counts: "dict[str, tuple[int]]") -> dict:
    stats = {}
    for sample, fn_list in infiles.items():
        stats[sample] = {'SampleID': sample}
        logger.info(f"Working on sample {sample} ...")

         # Process flagstats file:       
        flagstat_file = fn_list[0]
        logger.info(f"Processing flagstat file {flagstat_file} for sample {sample}.")
        flagstats = read_flagstats(flagstat_file)
        stats[sample].update(flagstats)

        # Process mosdepth summary file:
        summary_file = fn_list[1]
        logger.info(f"Processing mosdepth summary {summary_file} for sample {sample}.")
        depth = read_mosdepth(summary_file)
        stats[sample].update(depth)

        # Process mosdepth coverage file:
        coverage_file = fn_list[2]
        logger.info(f"Processing mosdepth coverage histogram {coverage_file} for sample {sample}.")
        coverage = read_coverage(coverage_file)
        stats[sample].update(coverage)

        # Add read counts:
        stats[sample]['TotalReads'] = counts[sample][0]
        stats[sample]['TrimmedReads'] = counts[sample][1]
        stats[sample]['%Trimmed'] = "N/A" if counts[sample][0] == "N/A" or counts[sample][0] == 0 \
                                        else f"{(counts[sample][1] / counts[sample][0]):.4f}"
        stats[sample]['%Filtered'] = "N/A" if counts[sample][1] == "N/A" or counts[sample][1] == 0 \
                                        else f"{(stats[sample]['FilteredReads'] / counts[sample][1]):.4f}"

    return stats


def write_output(stats: dict, outfile: Path, header: list):
    with open(outfile, 'w', newline='') as outhandle:
        writer = csv.DictWriter(outhandle, header, delimiter='\t')
        writer.writeheader()
        writer.writerows(stats.values())


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate YAML file for ERGA EAR.",
        epilog="Example: python generate_yaml.py template.yaml -f flagstat.txt -s summary.txt -c coverage.txt -o stats_summary.tsv",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "template",
        type=Path,
        help="YAML template file."
    )
    parser.add_argument(
        "-g",
        "--genomescope",
        type=Path,
        required=True,
        help="GenomeScope summary file"
    )
    parser.add_argument(
        "--smudgeplot",
        type=Path,
        required=False,
        help="Smudgeplot summary file"
    )
    parser.add_argument(
        "--merqury",
        type=Path,
        required=True,
        help="Merqury result folder"
    )
    parser.add_argument(
        "--hifi",
        type=Path,
        required=False,
        help="Mosdepth file for PacBio HiFi reads"
    )
    parser.add_argument(
        "--ul",
        type=Path,
        required=False,
        help="Mosdepth file for ONT UL reads"
    )
    parser.add_argument(
        "--ul",
        type=Path,
        required=False,
        help="Mosdepth file for HiC reads"
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=Path,
        default="input.yaml",
        help="YAML file for make_EAR.py"
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
    infiles, samples = get_infiles(args.flagstat, args.summary, args.coverage)
    try:
        check_infiles(infiles)
    except Exception as e:
        logger.error(e)
        sys.exit(1)
    try:
        counts = get_counts(samples, args.before, args.after)
    except Exception as e:
        logger.error(e)
        sys.exit(1)

    with open(args.template, 'r') as inhandle:
        data = yaml.safe_load(inhandle)
    
    data['']

    stats = read_stats(infiles, counts)
    write_output(stats, args.outfile, HEADER)


if __name__ == "__main__":
    sys.exit(main())
