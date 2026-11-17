"""Minimal tests for the run summary report: parsing and section selection."""

from pathlib import Path

from scripts.report_summary import load_stats, read_tsv, render_report

STATS = "mapped\t181\nmean coverage\t12\n% recovery >= 10x\t97\n"
BED = "start\tstop\tname\tmean coverage\t% recovery >= 10x\n1\t9\tamp_1\t5\t80\n"


def _write_inputs(root: Path, *, with_svg: bool = False, with_bed: bool = False) -> Path:
	stage = root / "sample" / ".bamdash"
	stage.mkdir(parents=True)
	(stage / "sample_plot_ref_bam_stats.tabular").write_text(STATS)
	if with_svg:
		(stage / "sample_plot_ref.svg").write_text("<svg></svg>")
	if with_bed:
		(stage / "sample_plot_ref_bed_data_1.tabular").write_text(BED)
	variants = root / "sample_variants.tsv"
	variants.write_text("CHROM\tPOS\nref\t42\n")
	return variants


def test_parses_stats_and_tsv(tmp_path: Path) -> None:
	stats_path = tmp_path / "stats.tabular"
	stats_path.write_text(STATS)
	stats, recovery_key = load_stats(stats_path)
	assert stats["mapped"] == "181"
	assert recovery_key == "% recovery >= 10x"

	tsv_path = tmp_path / "table.tsv"
	tsv_path.write_text(BED)
	headers, rows = read_tsv(tsv_path)
	assert headers[0] == "start"
	assert rows == [["1", "9", "amp_1", "5", "80"]]


def test_report_includes_and_excludes_sections(tmp_path: Path) -> None:
	output = tmp_path / "report_summary.html"
	variants = _write_inputs(tmp_path)
	stage = tmp_path / "sample" / ".bamdash"

	# Optional inputs absent -> their sections are excluded, nothing breaks.
	render_report({"sample": stage}, {"sample": variants}, {"sample": ["ref"]}, True, output)
	content = output.read_text()
	assert "181" in content and "97" in content and "42" in content
	assert "Coverage summary" in content and "Called variants" in content
	assert "data:image/svg+xml;base64," not in content
	assert "amp_1" not in content

	# Optional inputs present -> their sections are included.
	(stage / "sample_plot_ref.svg").write_text("<svg></svg>")
	(stage / "sample_plot_ref_bed_data_1.tabular").write_text(BED)
	render_report({"sample": stage}, {"sample": variants}, {"sample": ["ref"]}, True, output)
	content = output.read_text()
	assert "data:image/svg+xml;base64," in content
	assert "amp_1" in content
	assert "data:image/svg+xml;base64," in content
