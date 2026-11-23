#!/usr/bin/env python3
"""Generate length- and GC-matched, mutually non-overlapping background regions."""

from __future__ import annotations

import argparse
from bisect import bisect_left, bisect_right
import json
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pysam
from scipy.stats import ks_2samp, ttest_rel
from tqdm.auto import tqdm


PathLike = Union[str, Path]
DEFAULT_GENOME = Path("/home/minj/Genomes/ENCODE/hg38/hg38.fa")
DEFAULT_BLACKLIST = Path(
    "/data1/minj/footprint/iMab_inter_TF_model/0.rawdata/blacklist.bed"
)


def _read_bed(path: PathLike, min_columns: int = 3):
    records = []
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < min_columns:
                raise ValueError(f"{path}:{line_number}: expected at least {min_columns} columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid BED coordinates") from error
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: require 0 <= start < end")
            records.append(
                {
                    "chrom": fields[0],
                    "start": start,
                    "end": end,
                    "fields": fields,
                    "line_number": line_number,
                }
            )
    if not records:
        raise ValueError(f"No BED records found in {path}")
    return records


def _gc_fraction(sequence: str) -> float:
    if not sequence:
        raise ValueError("Cannot calculate GC for an empty sequence")
    sequence = sequence.upper()
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def _reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhd bn".replace(" ", ""))
    return sequence.translate(table)[::-1]


def _merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return merged


class ExclusionIndex:
    """Small dynamic interval index for disjoint BED intervals."""

    def __init__(self, intervals_by_chrom):
        self.starts = {}
        self.ends = {}
        for chrom, intervals in intervals_by_chrom.items():
            merged = _merge_intervals(intervals)
            self.starts[chrom] = [item[0] for item in merged]
            self.ends[chrom] = [item[1] for item in merged]

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        starts = self.starts.get(chrom, [])
        ends = self.ends.get(chrom, [])
        index = bisect_right(starts, end - 1) - 1
        return index >= 0 and ends[index] > start

    def add(self, chrom: str, start: int, end: int):
        if self.overlaps(chrom, start, end):
            raise ValueError(f"Cannot add overlapping interval {chrom}:{start}-{end}")
        starts = self.starts.setdefault(chrom, [])
        ends = self.ends.setdefault(chrom, [])
        index = bisect_left(starts, start)
        starts.insert(index, start)
        ends.insert(index, end)


def _validate_reference(records, fasta):
    references = set(fasta.references)
    lengths = dict(zip(fasta.references, fasta.lengths))
    for record in records:
        chrom = record["chrom"]
        if chrom not in references:
            raise ValueError(f"Chromosome {chrom!r} is absent from the reference FASTA")
        if record["end"] > lengths[chrom]:
            raise ValueError(
                f"Input interval {chrom}:{record['start']}-{record['end']} exceeds chromosome length"
            )
    return lengths


def _build_exclusions(input_records, chromosome_lengths, blacklist_bed, exclude_bed, min_distance):
    intervals = {}
    for record in input_records:
        chrom = record["chrom"]
        intervals.setdefault(chrom, []).append(
            (
                max(0, record["start"] - min_distance),
                min(chromosome_lengths[chrom], record["end"] + min_distance),
            )
        )

    for path in [blacklist_bed, exclude_bed]:
        if path is None:
            continue
        for record in _read_bed(path):
            chrom = record["chrom"]
            if chrom not in chromosome_lengths:
                continue
            start = max(0, record["start"])
            end = min(chromosome_lengths[chrom], record["end"])
            if end > start:
                intervals.setdefault(chrom, []).append((start, end))
    return ExclusionIndex(intervals)


def generate_gc_matched_background(
    input_bed: PathLike,
    output_prefix: PathLike,
    genome_fasta: PathLike = DEFAULT_GENOME,
    blacklist_bed: Optional[PathLike] = DEFAULT_BLACKLIST,
    exclude_bed: Optional[PathLike] = None,
    gc_tolerance: float = 0.02,
    max_attempts: int = 10_000,
    min_distance: int = 0,
    random_seed: int = 42,
    reject_ambiguous: bool = True,
    show_progress: bool = True,
):
    """Generate one background interval for every input BED interval.

    Guarantees, when the function returns successfully:
      * output and input counts are identical;
      * each output has the same length and chromosome as its paired input;
      * per-pair absolute GC difference is <= ``gc_tolerance``;
      * outputs do not overlap the input, blacklist, optional exclusions, or
        one another;
      * BED and FASTA outputs describe the same sequences.

    The input may have three or more BED columns. If column 6 is '+' or '-',
    that strand is retained; otherwise '+' is used.
    """
    input_bed = Path(input_bed).resolve()
    output_prefix = Path(output_prefix).resolve()
    genome_fasta = Path(genome_fasta).resolve()
    blacklist_bed = Path(blacklist_bed).resolve() if blacklist_bed else None
    exclude_bed = Path(exclude_bed).resolve() if exclude_bed else None

    if not 0 <= gc_tolerance <= 1:
        raise ValueError("gc_tolerance must be between 0 and 1")
    if max_attempts < 1 or min_distance < 0:
        raise ValueError("max_attempts must be positive and min_distance non-negative")
    if not genome_fasta.exists():
        raise FileNotFoundError(genome_fasta)
    if not Path(str(genome_fasta) + ".fai").exists():
        raise FileNotFoundError(f"FASTA index is missing: {genome_fasta}.fai")

    records = _read_bed(input_bed)
    fasta = pysam.FastaFile(str(genome_fasta))
    try:
        chromosome_lengths = _validate_reference(records, fasta)
        exclusions = _build_exclusions(
            records, chromosome_lengths, blacklist_bed, exclude_bed, min_distance
        )

        for index, record in enumerate(records):
            sequence = fasta.fetch(record["chrom"], record["start"], record["end"]).upper()
            record["index"] = index
            record["length"] = record["end"] - record["start"]
            record["input_sequence"] = sequence
            record["input_gc"] = _gc_fraction(sequence)

        # Long and GC-extreme intervals are placed first because they have fewer choices.
        placement_order = sorted(
            records,
            key=lambda item: (abs(item["input_gc"] - 0.42), item["length"]),
            reverse=True,
        )
        rng = np.random.default_rng(random_seed)
        results = [None] * len(records)

        progress = tqdm(
            placement_order,
            total=len(placement_order),
            desc="GC-matched background",
            unit="region",
            dynamic_ncols=True,
            disable=not show_progress,
        )
        maximum_attempts_seen = 0
        for record in progress:
            chrom = record["chrom"]
            length = record["length"]
            maximum_start = chromosome_lengths[chrom] - length
            if maximum_start < 0:
                raise RuntimeError(f"Interval length {length} does not fit on {chrom}")

            best = None
            for attempt in range(1, max_attempts + 1):
                start = int(rng.integers(0, maximum_start + 1))
                end = start + length
                if exclusions.overlaps(chrom, start, end):
                    continue
                sequence = fasta.fetch(chrom, start, end).upper()
                if reject_ambiguous and any(base not in "ACGT" for base in sequence):
                    continue
                gc = _gc_fraction(sequence)
                difference = abs(gc - record["input_gc"])
                if best is None or difference < best[0]:
                    best = (difference, start, end, gc, sequence, attempt)
                if difference <= gc_tolerance:
                    break

            if best is None or best[0] > gc_tolerance:
                best_text = "no valid location" if best is None else f"best GC difference={best[0]:.6f}"
                raise RuntimeError(
                    f"Could not match input line {record['line_number']} after {max_attempts} "
                    f"attempts ({best_text}). Increase gc_tolerance or max_attempts."
                )

            difference, start, end, gc, sequence, _best_attempt = best
            attempts_used = attempt
            exclusions.add(chrom, start, end)
            strand = record["fields"][5] if len(record["fields"]) >= 6 else "+"
            if strand not in {"+", "-"}:
                strand = "+"
            results[record["index"]] = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": f"bg_{record['index'] + 1}|source={chrom}:{record['start']}-{record['end']}",
                "strand": strand,
                "sequence": sequence,
                "input_gc": record["input_gc"],
                "background_gc": gc,
                "gc_difference": difference,
                "attempts": attempts_used,
            }
            if attempts_used > maximum_attempts_seen:
                maximum_attempts_seen = attempts_used
                progress.set_postfix(max_attempts=maximum_attempts_seen, refresh=False)
    finally:
        fasta.close()

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    bed_path = Path(str(output_prefix) + ".bed")
    fasta_path = Path(str(output_prefix) + ".fa")
    report_path = Path(str(output_prefix) + ".validation.json")

    bed_tmp = Path(str(bed_path) + ".tmp")
    fasta_tmp = Path(str(fasta_path) + ".tmp")
    with bed_tmp.open("w") as bed_handle, fasta_tmp.open("w") as fasta_handle:
        for result in results:
            bed_handle.write(
                f"{result['chrom']}\t{result['start']}\t{result['end']}\t"
                f"{result['name']}\t.\t{result['strand']}\n"
            )
            sequence = result["sequence"]
            if result["strand"] == "-":
                sequence = _reverse_complement(sequence)
            fasta_handle.write(f">{result['name']} strand={result['strand']}\n")
            for offset in range(0, len(sequence), 80):
                fasta_handle.write(sequence[offset : offset + 80] + "\n")
    bed_tmp.replace(bed_path)
    fasta_tmp.replace(fasta_path)

    input_gc = np.array([result["input_gc"] for result in results])
    background_gc = np.array([result["background_gc"] for result in results])
    differences = np.abs(input_gc - background_gc)
    input_only = ExclusionIndex(
        {
            chrom: [(record["start"], record["end"]) for record in records if record["chrom"] == chrom]
            for chrom in {record["chrom"] for record in records}
        }
    )
    input_overlap_count = sum(
        input_only.overlaps(result["chrom"], result["start"], result["end"])
        for result in results
    )
    background_self_overlap_count = 0
    for chrom in {result["chrom"] for result in results}:
        chrom_results = sorted(
            (result["start"], result["end"])
            for result in results
            if result["chrom"] == chrom
        )
        previous_end = -1
        for start, end in chrom_results:
            if start < previous_end:
                background_self_overlap_count += 1
            previous_end = max(previous_end, end)
    if input_overlap_count or background_self_overlap_count:
        raise AssertionError(
            f"Post-write interval validation failed: input overlaps={input_overlap_count}, "
            f"background self-overlaps={background_self_overlap_count}"
        )
    report = {
        "input_bed": str(input_bed),
        "genome_fasta": str(genome_fasta),
        "blacklist_bed": str(blacklist_bed) if blacklist_bed else None,
        "additional_exclude_bed": str(exclude_bed) if exclude_bed else None,
        "output_bed": str(bed_path),
        "output_fasta": str(fasta_path),
        "records": len(results),
        "same_chromosome_per_pair": True,
        "same_length_per_pair": True,
        "input_overlap_count": input_overlap_count,
        "background_self_overlap_count": background_self_overlap_count,
        "gc_tolerance": gc_tolerance,
        "maximum_absolute_gc_difference": float(differences.max()),
        "mean_absolute_gc_difference": float(differences.mean()),
        "input_gc_mean": float(input_gc.mean()),
        "background_gc_mean": float(background_gc.mean()),
        "paired_ttest_pvalue": float(ttest_rel(input_gc, background_gc).pvalue),
        "ks_test_pvalue": float(ks_2samp(input_gc, background_gc).pvalue),
        "random_seed": random_seed,
        "max_attempts_per_record": max_attempts,
        "maximum_attempts_used": int(max(result["attempts"] for result in results)),
        "min_distance_from_input": min_distance,
    }
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return report


def _optional_path(value):
    return None if value.lower() in {"none", "null", "no"} else value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_bed", help="Input BED (three or more columns)")
    parser.add_argument("output_prefix", help="Output prefix; .bed/.fa/.validation.json are added")
    parser.add_argument("--genome", default=str(DEFAULT_GENOME), help="Indexed reference FASTA")
    parser.add_argument(
        "--blacklist",
        type=_optional_path,
        default=str(DEFAULT_BLACKLIST),
        help="Blacklist BED, or 'none' (default: hg38 blacklist)",
    )
    parser.add_argument("--exclude-bed", type=_optional_path, default=None)
    parser.add_argument("--gc-tolerance", type=float, default=0.02)
    parser.add_argument("--max-attempts", type=int, default=10_000)
    parser.add_argument("--min-distance", type=int, default=0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--allow-ambiguous", action="store_true")
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable the progress bar (useful for batch logs)",
    )
    args = parser.parse_args()
    generate_gc_matched_background(
        input_bed=args.input_bed,
        output_prefix=args.output_prefix,
        genome_fasta=args.genome,
        blacklist_bed=args.blacklist,
        exclude_bed=args.exclude_bed,
        gc_tolerance=args.gc_tolerance,
        max_attempts=args.max_attempts,
        min_distance=args.min_distance,
        random_seed=args.seed,
        reject_ambiguous=not args.allow_ambiguous,
        show_progress=not args.no_progress,
    )


if __name__ == "__main__":
    main()
