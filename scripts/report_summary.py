#!/usr/bin/env python3
"""Prepare and render the self-contained run summary report."""

from __future__ import annotations

import argparse
import base64
import csv
import json
import re
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, StrictUndefined


def read_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
	"""Read a headered tab-separated file without changing its values."""
	with path.open(newline="", encoding="utf-8") as handle:
		rows = list(csv.reader(handle, delimiter="\t"))
	return (rows[0], rows[1:]) if rows else ([], [])


def load_stats(path: Path) -> tuple[dict[str, str], str]:
	"""Load a BAMdash stats TSV and locate its recovery column key."""
	with path.open(newline="", encoding="utf-8") as handle:
		stats = {row[0].strip(): row[1].strip() for row in csv.reader(handle, delimiter="\t") if len(row) >= 2}
	recovery_keys = [key for key in stats if key.startswith("% recovery >=")]
	if len(recovery_keys) != 1:
		raise ValueError(f"Expected exactly one recovery column in {path}, found {recovery_keys}")
	return stats, recovery_keys[0]


def _find_sidecar(stage_dir: Path, sample: str, reference: str, suffix: str, single: bool) -> Path | None:
	"""Find a BAMdash sidecar using its multi- or single-reference naming."""
	safe_ref = re.sub(r"[^A-Za-z0-9]+", "_", reference).strip("_")
	names = [f"{sample}_plot_{safe_ref}{suffix}"]
	if single:
		names.append(f"{sample}_plot{suffix}")
	return next((stage_dir / name for name in names if (stage_dir / name).is_file()), None)


def _amplicon_table(bed_path: Path, reference: str, recovery_key: str) -> dict:
	"""Select the reported amplicon columns from a BAMdash BED dump."""
	headers, rows = read_tsv(bed_path)
	selected = [col for col in ("name", "mean coverage", recovery_key) if col in headers]
	indices = [headers.index(col) for col in selected]
	return {
		"reference": reference,
		"headers": selected,
		"rows": [[row[i] for i in indices] for row in rows],
	}


def render_report(
		stage_dirs: dict[str, Path],
		variant_paths: dict[str, Path],
		sample_references: dict[str, list[str]],
		include_svg: bool,
		output_path: Path,
		config_yaml: str | None = None,
		sample_sheet: str | None = None,
) -> None:
	"""Read sample outputs and write a self-contained HTML report."""
	summary: list[list[str]] = []
	samples: list[dict] = []
	recovery_key: str | None = None

	for sample, references in sample_references.items():
		stage = stage_dirs[sample]
		single = len(references) == 1
		plots: list[dict] = []
		amplicons: list[dict] = []

		for reference in references:
			stats_path = _find_sidecar(stage, sample, reference, "_bam_stats.tabular", single)
			if stats_path is None:
				raise FileNotFoundError(f"Missing BAMdash stats for sample {sample}, reference {reference}")
			stats, recovery_key = load_stats(stats_path)
			summary.append([sample, reference, stats["mapped"], stats["mean coverage"], stats[recovery_key]])

			svg_path = _find_sidecar(stage, sample, reference, ".svg", single) if include_svg else None
			if svg_path:
				plots.append({
					"reference": reference,
					"image": base64.b64encode(svg_path.read_bytes()).decode("ascii"),
				})

			bed_path = _find_sidecar(stage, sample, reference, "_bed_data_1.tabular", single)
			if bed_path:
				amplicons.append(_amplicon_table(bed_path, reference, recovery_key))

		variant_headers, variant_rows = read_tsv(variant_paths[sample])
		samples.append({
			"name": sample,
			"plots": plots,
			"amplicons": amplicons,
			"variants": {"headers": variant_headers, "rows": variant_rows},
		})

	if recovery_key is None:
		raise ValueError("No BAMdash stats found for any sample")

	env = Environment(
		loader=FileSystemLoader(Path(__file__).resolve().parent / "report_summary"),
		autoescape=True,
		undefined=StrictUndefined,
	)
	document = env.get_template("report.html.j2").render(
		summary=summary,
		recovery_key=recovery_key,
		samples=samples,
		config_yaml=config_yaml.strip() if config_yaml else None,
		sample_sheet=sample_sheet.strip() if sample_sheet else None,
	)
	output_path.write_text(document, encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("--bamdash-dirs", type=Path, required=True)
	parser.add_argument("--variants", type=Path, required=True)
	parser.add_argument("--references", type=Path, required=True)
	parser.add_argument("--config", type=Path, default=None)
	parser.add_argument("--sample-sheet", type=Path, default=None)
	parser.add_argument("--svg", action=argparse.BooleanOptionalAction, default=True)
	parser.add_argument("-o", "--output", type=Path, required=True)
	return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
	args = parse_args(argv)
	bamdash_paths = [Path(line) for line in args.bamdash_dirs.read_text(encoding="utf-8").splitlines()]
	variant_paths = [Path(line) for line in args.variants.read_text(encoding="utf-8").splitlines()]
	references = json.loads(args.references.read_text(encoding="utf-8"))
	config_yaml = args.config.read_text(encoding="utf-8") if args.config else None
	sample_sheet = args.sample_sheet.read_text(encoding="utf-8") if args.sample_sheet else None
	samples = list(references)

	render_report(
		dict(zip(samples, bamdash_paths, strict=True)),
		dict(zip(samples, variant_paths, strict=True)),
		references,
		args.svg,
		args.output,
		config_yaml=config_yaml,
		sample_sheet=sample_sheet,
	)


if __name__ == "__main__":
	main()
