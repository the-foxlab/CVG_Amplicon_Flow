"""Utility functions shared by the Snakemake workflow.

These helpers keep parsing and hashing logic out of the `Snakefile` so the
workflow rules stay focused on rule definitions.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Mapping

import pandas as pd
from Bio import SeqIO


def hashfile(file_names: list[str]) -> str:
    """Compute a SHA256 digest across one or more files.

    Files are read in chunks to avoid loading large references into memory.
    """

    buffer_size = 65536
    sha256 = hashlib.sha256()

    for file_name in file_names:
        with open(file_name, "rb") as handle:
            while True:
                data = handle.read(buffer_size)
                if not data:
                    break
                sha256.update(data)

    return sha256.hexdigest()


def is_missing(value: Any) -> bool:
    """Return True if a value should be treated as missing."""

    text = str(value).strip()
    return not text or text.lower() == "nan"


def read_fasta_names(fasta_path: str) -> list[str]:
    """Return unique FASTA record IDs from a FASTA file."""

    names: list[str] = []
    seen: set[str] = set()
    with open(fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            name = record.id
            if name in seen:
                raise ValueError(f"Duplicate FASTA record name '{name}' in {fasta_path}")
            seen.add(name)
            names.append(name)
    if not names:
        raise ValueError(f"No FASTA records found in {fasta_path}")
    return names


def read_bed_names(bed_path: str) -> list[str]:
    """Read unique reference names from the first BED column."""

    names: list[str] = []
    seen: set[str] = set()
    with open(bed_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0].strip()
            if not name:
                continue
            if name not in seen:
                seen.add(name)
                names.append(name)
    return names


def load_samples_from_sheet(sample_sheet: str) -> dict[str, dict[str, Any]]:
    """Load and validate the pipeline sample sheet into an internal map."""

    samples: dict[str, dict[str, Any]] = {}
    table = pd.read_csv(sample_sheet, sep=",", engine="python")

    required_columns = {"sample", "reference_fasta", "r1", "r2"}
    if not required_columns.issubset(set(table.columns)):
        raise ValueError("sample_sheet must contain columns: sample, reference_fasta, r1, r2")

    has_bed_column = "bed" in set(table.columns)

    for row in table.itertuples(index=False):
        sample = str(row.sample).strip()
        if is_missing(sample):
            continue

        reference_fasta = str(row.reference_fasta).strip()
        r1 = str(row.r1).strip()
        r2 = str(row.r2).strip()
        bed = str(row.bed).strip() if has_bed_column else ""

        if is_missing(reference_fasta) or is_missing(r1) or is_missing(r2):
            raise ValueError(f"Sample '{sample}' has missing reference_fasta, r1 or r2 in sample_sheet")

        if sample in samples:
            raise ValueError(f"Duplicate sample in sample_sheet: '{sample}'")

        reference_names = read_fasta_names(reference_fasta)
        bed_names: list[str] = []
        if not is_missing(bed):
            bed_names = read_bed_names(bed)
            unknown_bed_names = sorted(set(bed_names) - set(reference_names))
            if unknown_bed_names:
                raise ValueError(
                    f"BED file for sample '{sample}' contains references not present in FASTA: {unknown_bed_names}"
                )

        effective_references = bed_names if bed_names else reference_names
        if not effective_references:
            raise ValueError(f"Sample '{sample}' has no usable references after FASTA/BED parsing")

        samples[sample] = {
            "reference_fasta": reference_fasta,
            "bed": "" if is_missing(bed) else bed,
            "r1": r1,
            "r2": r2,
            "references": effective_references,
            "reference_hash": hashfile([reference_fasta]) if is_missing(bed) else hashfile([reference_fasta, bed]),
        }

    if not samples:
        raise ValueError("sample_sheet was provided but no samples were parsed")
    return samples


def load_samples(cfg: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    """Load sample map from config values."""

    sample_sheet = cfg.get("sample_sheet")
    if not sample_sheet:
        raise ValueError("sample_sheet is required")
    return load_samples_from_sheet(sample_sheet)


def append_res_dir(res_dir: str, path_fragment: str) -> str:
    """Join a path fragment under the workflow result directory."""

    return str(Path(res_dir) / path_fragment)


def build_bwa_files(sample_entry: Mapping[str, Any], res_dir: str) -> dict[str, str]:
    """Build all reference index paths needed by BWA and samtools."""

    hash_value = str(sample_entry["reference_hash"])
    ref_type = "masked" if sample_entry["bed"] else "unmasked"
    ref = append_res_dir(res_dir, f"references/{hash_value}/{hash_value}.{ref_type}.fasta")
    return {
        "ref": ref,
        "ref_fai": f"{ref}.fai",
        "ref_amb": f"{ref}.amb",
        "ref_ann": f"{ref}.ann",
        "ref_bwt": f"{ref}.bwt",
        "ref_pac": f"{ref}.pac",
        "ref_sa": f"{ref}.sa",
    }


def masked_reference_fasta(sample_entry: Mapping[str, Any], res_dir: str) -> str:
    """Return the masked/unmasked reference path for a sample."""

    hash_value = str(sample_entry["reference_hash"])
    ref_type = "masked" if sample_entry["bed"] else "unmasked"
    return append_res_dir(res_dir, f"references/{hash_value}/{hash_value}.{ref_type}.fasta")