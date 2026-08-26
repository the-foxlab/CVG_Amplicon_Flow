#!/usr/bin/env python3

import argparse
import logging
from collections import defaultdict
from pathlib import Path
from typing import NamedTuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


logger = logging.getLogger(__name__)


class Primer(NamedTuple):
    reference: str
    start: int
    end: int


def read_primers(fn_bed: Path) -> dict[str, list[Primer]]:
    primers = defaultdict(list)

    with fn_bed.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 3:
                raise ValueError(
                    f"Invalid BED line {line_number}: expected at least 3 columns."
                )

            try:
                primer = Primer(
                    reference=fields[0],
                    start=int(fields[1]),
                    end=int(fields[2]),
                )
            except ValueError as exc:
                raise ValueError(
                    f"Invalid BED coordinates on line {line_number}: {line}"
                ) from exc

            if primer.start < 0 or primer.end <= primer.start:
                raise ValueError(
                    f"Invalid BED interval on line {line_number}: "
                    f"{primer.start}-{primer.end}"
                )

            primers[primer.reference].append(primer)

    logger.info(
        "Read primers for %d reference(s) from %s",
        len(primers),
        fn_bed,
    )

    return dict(primers)


def mask_fasta_beginnings(
    fasta_obj: SeqRecord,
    primers: list[Primer],
) -> SeqRecord:
    if not primers:
        return fasta_obj

    first_primer = min(primers, key=lambda primer: primer.start)
    last_primer = max(primers, key=lambda primer: primer.end)

    ref_len =  len(fasta_obj.seq)

    if first_primer.end > ref_len:
        raise ValueError(
            f"Primer end {first_primer.end} exceeds sequence length "
            f"{len(fasta_obj.seq)} for '{fasta_obj.id}'."
        )

    if last_primer.start > ref_len:
        raise ValueError(
            f"Primer start {last_primer.start} exceeds sequence length "
            f"{len(fasta_obj.seq)} for '{fasta_obj.id}'."
        )

    logger.info(
        "Masking %s: 0-%d and %d-%d",
        fasta_obj.id,
        first_primer.end,
        last_primer.start,
        len(fasta_obj.seq)
    )

    fasta_obj.seq = Seq(
        "N" * first_primer.end
        + str(fasta_obj.seq[first_primer.end:last_primer.start])
        + "N" * (len(fasta_obj.seq) - last_primer.start)
    )

    return fasta_obj


def prepare_fastas(
    fn_fasta: Path,
    fn_bed: Path,
):
    primers = read_primers(fn_bed)

    for fasta_obj in SeqIO.parse(fn_fasta, "fasta"):
        reference_primers = primers.get(fasta_obj.id)

        if reference_primers is None:
            logger.warning(
                "No primers found for reference '%s'. Sequence unchanged.",
                fasta_obj.id,
            )
        else:
            fasta_obj = mask_fasta_beginnings(
                fasta_obj,
                reference_primers,
            )

        yield fasta_obj


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Mask FASTA sequence ends according to the outermost primers "
            "defined in a BED file."
        )
    )

    parser.add_argument(
        "-f",
        "--fasta",
        type=Path,
        required=True,
        help="Input FASTA file.",
    )

    parser.add_argument(
        "-b",
        "--bed",
        type=Path,
        required=True,
        help="BED file containing primer coordinates.",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output FASTA file.",
    )

    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    logger.info("Processing FASTA file %s", args.fasta)

    n_written = SeqIO.write(
        prepare_fastas(args.fasta, args.bed),
        args.output,
        "fasta",
    )

    logger.info(
        "Wrote %d sequence(s) to %s",
        n_written,
        args.output,
    )


if __name__ == "__main__":
    main()
