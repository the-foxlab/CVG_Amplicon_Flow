#!/usr/bin/env python3
"""Convert a primer BED file into an amplicon BED file.

The input BED lists primer binding sites in alternating order: a ``_LEFT``
primer followed by its matching ``_RIGHT`` primer, one amplicon per pair:

    HSV_1_UL23	0	18	UL23_0_LEFT	1	+
    HSV_1_UL23	485	506	UL23_0_RIGHT	1	-

Each pair is merged into a single amplicon interval spanning both primers:

    HSV_1_UL23	0	506	UL23_0	1	+

Amplicons are emitted with score 1 and strand ``+``.

The output is suitable as an additional track for the bamdash plot.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import NamedTuple

logger = logging.getLogger(__name__)

LEFT_SUFFIX = "_LEFT"
RIGHT_SUFFIX = "_RIGHT"


class BedEntry(NamedTuple):
    reference: str
    start: int
    end: int
    name: str
    score: int
    strand: str


def read_primer_bed(bed_path: Path) -> list[BedEntry]:
    """Read a 6-column primer BED file."""

    entries: list[BedEntry] = []

    with bed_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 3:
                raise ValueError(
                    f"Invalid BED line {line_number}: expected at least 3 columns."
                )

            if len(fields) < 6:
                raise ValueError(
                    f"Invalid BED line {line_number}: primer BED requires 6 columns "
                    f"(reference, start, end, name, score, strand)."
                )

            try:
                entry = BedEntry(
                    reference=fields[0],
                    start=int(fields[1]),
                    end=int(fields[2]),
                    name=fields[3],
                    score=int(fields[4]),
                    strand=fields[5],
                )
            except ValueError as exc:
                raise ValueError(
                    f"Invalid BED fields on line {line_number}: {line}"
                ) from exc

            if entry.start < 0 or entry.end <= entry.start:
                raise ValueError(
                    f"Invalid BED interval on line {line_number}: "
                    f"{entry.start}-{entry.end}"
                )

            entries.append(entry)

    if not entries:
        raise ValueError(f"No primer entries found in {bed_path}")

    return entries


def primer_to_amplicon_name(name: str, line_number: int) -> str:
    """Strip the ``_LEFT``/``_RIGHT`` suffix from a primer name."""

    if name.endswith(LEFT_SUFFIX):
        return name[: -len(LEFT_SUFFIX)]
    if name.endswith(RIGHT_SUFFIX):
        return name[: -len(RIGHT_SUFFIX)]

    raise ValueError(
        f"Primer name '{name}' on line {line_number} does not end with "
        f"'{LEFT_SUFFIX}' or '{RIGHT_SUFFIX}'."
    )


def build_amplicons(entries: list[BedEntry]) -> list[BedEntry]:
    """Merge consecutive LEFT/RIGHT primer pairs into amplicon intervals.

    Each amplicon is a LEFT primer immediately followed by its matching RIGHT
    primer. The amplicon spans from the leftmost primer start to the rightmost
    primer end and is emitted with score (defines pool) 1 and strand ``+``.
    """

    if len(entries) % 2 != 0:
        raise ValueError(
            "Primer BED must contain an even number of entries "
            f"(LEFT/RIGHT pairs), found {len(entries)}."
        )

    amplicons: list[BedEntry] = []

    for index in range(0, len(entries), 2):
        left = entries[index]
        right = entries[index + 1]
        line_number = index + 1

        left_name = primer_to_amplicon_name(left.name, line_number)
        right_name = primer_to_amplicon_name(right.name, line_number + 1)

        if right_name != left_name:
            raise ValueError(
                f"Mismatched primer pair on lines {line_number}-{line_number + 1}: "
                f"'{left.name}' and '{right.name}'."
            )

        if left.strand not in ("+", "-") or right.strand not in ("+", "-"):
            raise ValueError(
                f"Invalid strand on primer pair '{left.name}'/'{right.name}': "
                f"'{left.strand}'/'{right.strand}'."
            )

        amplicon = BedEntry(
            reference=left.reference,
            start=min(left.start, right.start),
            end=max(left.end, right.end),
            name=left_name,
            score=1,
            strand="+",
        )
        amplicons.append(amplicon)

    return amplicons


def write_amplicon_bed(amplicons: list[BedEntry], output_path: Path) -> None:
    """Write amplicon entries as a 6-column BED file."""

    with output_path.open("w") as handle:
        for amplicon in amplicons:
            handle.write(
                f"{amplicon.reference}\t{amplicon.start}\t{amplicon.end}\t"
                f"{amplicon.name}\t{amplicon.score}\t{amplicon.strand}\n"
            )


def convert_primer_bed(bed_path: Path, output_path: Path) -> list[BedEntry]:
    """Convert a primer BED file to an amplicon BED file."""

    entries = read_primer_bed(bed_path)
    amplicons = build_amplicons(entries)
    write_amplicon_bed(amplicons, output_path)

    logger.info(
        "Converted %d primer(s) into %d amplicon(s): %s -> %s",
        len(entries),
        len(amplicons),
        bed_path,
        output_path,
    )
    return amplicons


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a primer BED file into an amplicon BED file."
    )
    parser.add_argument("primer_bed", type=Path, help="Input primer BED file.")
    parser.add_argument("output_bed", type=Path, help="Output amplicon BED file.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    args = parse_args(argv)
    convert_primer_bed(args.primer_bed, args.output_bed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
