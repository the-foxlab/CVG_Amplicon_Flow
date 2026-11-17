"""Prepare and render the self-contained run summary report.

The report is built from BAMdash's dumped sidecars plus the sample variants
TSV. All data is collected first, then rendered through a Jinja template so
the generated HTML stays self-contained (CSS/JS inlined, SVGs embedded as
data images).
"""

from __future__ import annotations

import base64
import csv
import re
from pathlib import Path

from jinja2 import (
	Environment,
	FileSystemLoader,
	StrictUndefined,
	Template,
	select_autoescape,
)


def sanitize_reference(reference: str) -> str:
	"""Match BAMdash's reference identifier in per-reference filenames."""
	return re.sub(r"[^A-Za-z0-9]+", "_", reference).strip("_")


def read_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
	"""Read a headered tab-separated file without changing its values."""
	with path.open(newline="", encoding="utf-8") as handle:
		rows = list(csv.reader(handle, delimiter="\t"))
	return (rows[0], rows[1:]) if rows else ([], [])


def find_reference_file(
		stage_dir: Path, sample: str, reference: str, suffix: str, single: bool
) -> Path | None:
	"""Find a BAMdash sidecar using its single- or multi-reference naming."""
	stem = f"{sample}_plot_{sanitize_reference(reference)}"
	candidate = stage_dir / f"{stem}{suffix}"
	if candidate.is_file():
		return candidate
	if single:
		candidate = stage_dir / f"{sample}_plot{suffix}"
		if candidate.is_file():
			return candidate
	return None


def load_stats(path: Path) -> tuple[dict[str, str], str]:
	"""Load a stats TSV and locate its recovery column key."""
	with path.open(newline="", encoding="utf-8") as handle:
		stats = {row[0].strip(): row[1].strip() for row in csv.reader(handle, delimiter="\t") if len(row) >= 2}
	recovery_keys = [key for key in stats if re.fullmatch(r"% recovery >= \S+x", key)]
	if len(recovery_keys) != 1:
		raise ValueError(f"Expected exactly one recovery column in {path}, found {recovery_keys}")
	return stats, recovery_keys[0]


def _check_reference_names(sample: str, references: list[str]) -> None:
	"""Reject sanitized names that collide or are empty for multi-reference samples."""
	safe_names = [sanitize_reference(reference) for reference in references]
	if len(references) > 1 and (not all(safe_names) or len(set(safe_names)) != len(safe_names)):
		raise ValueError(f"BAMdash reference filename collision or empty identifier for sample {sample}")


def _unify_recovery_key(current: str | None, candidate: str) -> str:
	"""Adopt the first recovery key and require all others to match it."""
	if current is None:
		return candidate
	if current != candidate:
		raise ValueError(f"Inconsistent recovery thresholds across samples: {current} vs {candidate}")
	return current


def _sidecar_finder(stage_dir: Path, sample: str, reference: str, single: bool):
	"""Return a lookup function for one sample/reference's sidecar files."""
	def find(suffix: str) -> Path | None:
		return find_reference_file(stage_dir, sample, reference, suffix, single)
	return find


def _amplicon_table(bed_path: Path, reference: str, recovery_key: str) -> dict:
	"""Select the reported amplicon columns from a BED-data dump."""
	headers, rows = read_tsv(bed_path)
	selected = [name for name in ["name", "mean coverage", recovery_key] if name in headers]
	indices = [headers.index(name) for name in selected]
	return {
		"reference": reference,
		"headers": selected,
		"rows": [[row[index] if index < len(row) else "" for index in indices] for row in rows],
	}


def _template() -> Template:
	"""Load the autoescaping HTML template from the report assets."""
	environment = Environment(
		loader=FileSystemLoader(Path(__file__).resolve().parent / "report_summary"),
		autoescape=select_autoescape(default=True),
		undefined=StrictUndefined,
	)
	return environment.get_template("report.html.j2")


def render_report(
		stage_dirs: dict[str, Path],
		variant_paths: dict[str, Path],
		sample_references: dict[str, list[str]],
		include_svg: bool,
		output_path: Path,
) -> None:
	"""Read sample outputs and write a portable, self-contained HTML report."""
	summary: list[list[str]] = []
	coverage: list[dict] = []
	variants: list[dict] = []
	amplicons: list[dict] = []
	recovery_key: str | None = None

	for sample, references in sample_references.items():
		_check_reference_names(sample, references)
		stage = stage_dirs[sample]

		variant_headers, variant_rows = read_tsv(variant_paths[sample])
		variants.append({"sample": sample, "headers": variant_headers, "rows": variant_rows})

		plots: list[dict] = []
		amplicon_tables: list[dict] = []
		for reference in references:
			sidecar = _sidecar_finder(stage, sample, reference, len(references) == 1)

			stats_path = sidecar("_bam_stats.tabular")
			if stats_path is None:
				raise FileNotFoundError(f"Missing BAMdash stats for sample {sample}, reference {reference}")
			stats, reference_key = load_stats(stats_path)
			recovery_key = _unify_recovery_key(recovery_key, reference_key)
			if not {"mapped", "mean coverage"}.issubset(stats):
				raise ValueError(f"BAMdash stats missing required metrics in {stats_path}")
			summary.append(
				[sample, reference, stats["mapped"], stats["mean coverage"], stats[recovery_key]]
			)

			svg_path = sidecar(".svg") if include_svg else None
			if svg_path:
				plots.append({
					"reference": reference,
					"image": base64.b64encode(svg_path.read_bytes()).decode("ascii"),
				})

			bed_path = sidecar("_bed_data_1.tabular")
			if bed_path:
				amplicon_tables.append(_amplicon_table(bed_path, reference, recovery_key))

		if plots:
			coverage.append({"sample": sample, "plots": plots})
		if amplicon_tables:
			amplicons.append({"sample": sample, "tables": amplicon_tables})

	if recovery_key is None:
		raise ValueError("No BAMdash stats found for any sample")

	# One ordered sample list drives the single global selector; every sample
	# appears in the variants section, so it is a safe fallback for ordering.
	samples = list(dict.fromkeys(
		[item["sample"] for item in coverage]
		+ [item["sample"] for item in variants]
		+ [item["sample"] for item in amplicons]
	))

	document = _template().render(
		summary=summary,
		recovery_key=recovery_key,
		samples=samples,
		coverage=coverage,
		variants=variants,
		amplicons=amplicons,
	)
	output_path.write_text(document, encoding="utf-8")


