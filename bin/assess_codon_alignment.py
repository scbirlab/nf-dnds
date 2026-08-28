#!/usr/bin/env python3

from __future__ import annotations

import argparse
from collections import Counter
import csv
import textwrap
from pathlib import Path


VALID_DNA = set("ACGT")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Assess a codon alignment for suitability for phylogenetic "
            "dN/dS analysis."
        )
    )

    parser.add_argument("--alignment", required=True)
    parser.add_argument("--taxon-id", required=True)
    parser.add_argument("--orthogroup-id", required=True)

    parser.add_argument("--assessment", default="assessment.tsv")
    parser.add_argument("--groups", default="sequence_groups.tsv")
    parser.add_argument("--analysable", default="analysable.fna")

    parser.add_argument("--min-sequences", type=int, default=3)
    parser.add_argument("--min-unique-sequences", type=int, default=2)
    parser.add_argument("--min-codons", type=int, default=20)
    parser.add_argument("--min-variable-codons", type=int, default=1)
    parser.add_argument("--max-ambiguous-fraction", type=float, default=0.10)

    return parser.parse_args()


def read_fasta(path: Path):
    records = []
    seen_ids = set()

    current_id = None
    current_seq = []

    def flush():
        if current_id is None:
            return

        seq = "".join(current_seq).upper().replace("U", "T")

        if current_id in seen_ids:
            raise ValueError(f"Duplicate FASTA identifier: {current_id}")

        seen_ids.add(current_id)
        records.append((current_id, seq))

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                flush()

                current_id = line[1:].split()[0]

                if not current_id:
                    raise ValueError("Empty FASTA identifier")

                current_seq = []

            else:
                if current_id is None:
                    raise ValueError(
                        "Sequence data encountered before first FASTA header"
                    )

                current_seq.append("".join(line.split()))

    flush()

    if not records:
        raise ValueError("Alignment contains no sequences")

    return records


def validate_alignment(records):
    lengths = {len(seq) for _, seq in records}

    if len(lengths) != 1:
        raise ValueError(
            "Sequences have unequal lengths: "
            + ", ".join(map(str, sorted(lengths)))
        )

    length = next(iter(lengths))

    if length == 0:
        raise ValueError("Alignment has zero length")

    if length % 3 != 0:
        raise ValueError(
            f"Alignment length ({length}) is not divisible by three"
        )

    return length


def assess(records):
    alignment_length = validate_alignment(records)

    n_sequences = len(records)
    n_codons = alignment_length // 3

    sequences = [seq for _, seq in records]

    # Exact nucleotide haplotypes.
    sequence_groups = {}

    for name, seq in records:
        sequence_groups.setdefault(seq, []).append(name)

    n_unique_sequences = len(sequence_groups)
    n_duplicate_sequences = n_sequences - n_unique_sequences

    n_callable_codons = 0
    n_variable_codons = 0
    n_variable_nt_sites = 0
    n_parsimony_informative_codons = 0
    n_ambiguous_codon_cells = 0

    for start in range(0, alignment_length, 3):

        codons = [
            seq[start:start + 3]
            for seq in sequences
        ]

        ambiguous = [
            any(base not in VALID_DNA for base in codon)
            for codon in codons
        ]

        n_ambiguous_codon_cells += sum(ambiguous)

        # Do not count partially missing codon columns as evidence
        # for evolutionary variation.
        if any(ambiguous):
            continue

        n_callable_codons += 1

        codon_counts = Counter(codons)

        if len(codon_counts) > 1:
            n_variable_codons += 1

            # A codon column is parsimony informative if at least
            # two different states occur at least twice each.
            if sum(count >= 2 for count in codon_counts.values()) >= 2:
                n_parsimony_informative_codons += 1

        for offset in range(3):
            bases = {
                codon[offset]
                for codon in codons
            }

            if len(bases) > 1:
                n_variable_nt_sites += 1

    total_codon_cells = n_sequences * n_codons

    ambiguous_fraction = (
        n_ambiguous_codon_cells / total_codon_cells
        if total_codon_cells
        else 0.0
    )

    callable_fraction = (
        n_callable_codons / n_codons
        if n_codons
        else 0.0
    )

    return {
        "alignment_length_nt": alignment_length,
        "n_sequences": n_sequences,
        "n_unique_sequences": n_unique_sequences,
        "n_duplicate_sequences": n_duplicate_sequences,
        "n_codons": n_codons,
        "n_callable_codons": n_callable_codons,
        "callable_fraction": callable_fraction,
        "n_variable_codons": n_variable_codons,
        "n_variable_nt_sites": n_variable_nt_sites,
        "n_parsimony_informative_codons": n_parsimony_informative_codons,
        "n_ambiguous_codon_cells": n_ambiguous_codon_cells,
        "ambiguous_fraction": ambiguous_fraction,
        "sequence_groups": sequence_groups,
    }


def classify(stats, args):
    """
    Return one primary status plus all applicable QC flags.

    Priority is deliberate:
      structural/sample-quality limitations
      -> invariant
      -> low variation
      -> analysable
    """

    flags = []

    if stats["n_sequences"] < args.min_sequences:
        flags.append("TOO_FEW_SEQUENCES")

    if stats["n_codons"] < args.min_codons:
        flags.append("TOO_SHORT")

    if stats["ambiguous_fraction"] > args.max_ambiguous_fraction:
        flags.append("TOO_MUCH_MISSINGNESS")

    if stats["n_variable_codons"] == 0:
        flags.append("INVARIANT")

    elif (
        stats["n_unique_sequences"] < args.min_unique_sequences
        or stats["n_variable_codons"] < args.min_variable_codons
    ):
        flags.append("LOW_VARIATION")

    precedence = [
        "TOO_FEW_SEQUENCES",
        "TOO_SHORT",
        "TOO_MUCH_MISSINGNESS",
        "INVARIANT",
        "LOW_VARIATION",
    ]

    status = next(
        (flag for flag in precedence if flag in flags),
        "ANALYSABLE",
    )

    reasons = {
        "TOO_FEW_SEQUENCES": (
            f"{stats['n_sequences']} sequences; "
            f"minimum is {args.min_sequences}"
        ),
        "TOO_SHORT": (
            f"{stats['n_codons']} codons; "
            f"minimum is {args.min_codons}"
        ),
        "TOO_MUCH_MISSINGNESS": (
            f"ambiguous codon-cell fraction "
            f"{stats['ambiguous_fraction']:.4f}; "
            f"maximum is {args.max_ambiguous_fraction}"
        ),
        "INVARIANT": (
            "no variable fully callable codon positions"
        ),
        "LOW_VARIATION": (
            f"{stats['n_unique_sequences']} unique sequences and "
            f"{stats['n_variable_codons']} variable codons"
        ),
        "ANALYSABLE": (
            "alignment passes configured analysis thresholds"
        ),
    }

    return status, flags, reasons[status]


def write_assessment(path, taxon_id, orthogroup_id, status, flags, reason, stats):
    columns = [
        "taxon_id",
        "orthogroup_id",
        "status",
        "qc_flags",
        "reason",
        "n_sequences",
        "n_unique_sequences",
        "n_duplicate_sequences",
        "alignment_length_nt",
        "n_codons",
        "n_callable_codons",
        "callable_fraction",
        "n_variable_codons",
        "n_variable_nt_sites",
        "n_parsimony_informative_codons",
        "n_ambiguous_codon_cells",
        "ambiguous_fraction",
    ]

    row = {
        "taxon_id": taxon_id,
        "orthogroup_id": orthogroup_id,
        "status": status,
        "qc_flags": ";".join(flags),
        "reason": reason,
        **{
            key: value
            for key, value in stats.items()
            if key != "sequence_groups"
        },
    }

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(row)


def write_groups(path, groups):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")

        writer.writerow([
            "representative",
            "n_copies",
            "members",
        ])

        for members in groups.values():
            representative = members[0]

            writer.writerow([
                representative,
                len(members),
                ",".join(members),
            ])


def write_deduplicated_fasta(path, groups):
    """
    Write one representative sequence for each exact nucleotide
    haplotype.

    Exact duplicates contribute zero-length branches and no
    substitution information, so they need not enter FastTree/HyPhy.
    """

    with open(path, "w") as handle:
        for sequence, members in groups.items():
            representative = members[0]

            handle.write(f">{representative}\n")

            for line in textwrap.wrap(sequence, width=80):
                handle.write(f"{line}\n")


def main():
    args = parse_args()

    if args.min_sequences < 1:
        raise ValueError("--min-sequences must be >= 1")

    if args.min_unique_sequences < 1:
        raise ValueError("--min-unique-sequences must be >= 1")

    if args.min_codons < 1:
        raise ValueError("--min-codons must be >= 1")

    if args.min_variable_codons < 1:
        raise ValueError("--min-variable-codons must be >= 1")

    if not 0 <= args.max_ambiguous_fraction <= 1:
        raise ValueError("--max-ambiguous-fraction must be between 0 and 1")

    records = read_fasta(Path(args.alignment))
    stats = assess(records)

    status, flags, reason = classify(
        stats,
        args,
    )

    write_assessment(
        args.assessment,
        args.taxon_id,
        args.orthogroup_id,
        status,
        flags,
        reason,
        stats,
    )

    write_groups(
        args.groups,
        stats["sequence_groups"],
    )

    # The existence of this file is the Nextflow routing mechanism.
    # Non-analysable alignments simply do not emit this channel item.
    if status == "ANALYSABLE":
        write_deduplicated_fasta(
            args.analysable,
            stats["sequence_groups"],
        )

    print(
        f"{args.taxon_id}:{args.orthogroup_id}\t"
        f"{status}\t"
        f"n={stats['n_sequences']}\t"
        f"unique={stats['n_unique_sequences']}\t"
        f"codons={stats['n_codons']}\t"
        f"variable_codons={stats['n_variable_codons']}"
    )


if __name__ == "__main__":
    main()
