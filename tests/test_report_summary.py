"""Minimal tests for the run summary report: parsing and section selection."""

import json
from pathlib import Path

from scripts.report_summary import load_stats, main, read_tsv, render_report

STATS = "mapped\t181\nmean coverage\t12\n% recovery >= 10x\t97\n"
BED = "start\tstop\tname\tmean coverage\t% recovery >= 10x\n1\t9\tamp_1\t5\t80\n"


def _write_sample(root: Path, name: str = "sample") -> tuple[Path, Path]:
	stage = root / name / ".bamdash"
	stage.mkdir(parents=True)
	(stage / f"{name}_plot_ref_bam_stats.tabular").write_text(STATS)
	variants = root / f"{name}_variants.tsv"
	variants.write_text("CHROM\tPOS\nref\t42\n")
	return stage, variants


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
	stage, variants = _write_sample(tmp_path)

	render_report({"sample": stage}, {"sample": variants}, {"sample": ["ref"]}, True, output)
	content = output.read_text()
	assert "181" in content and "97" in content and "42" in content
	assert "Illumina sequencing · run report · Generated at: " in content
	assert "Coverage summary" in content and "Called variants" in content
	assert "select-sample" not in content
	assert "data:image/svg+xml;base64," not in content
	assert "amp_1" not in content
	assert "Run configuration" not in content
	assert "Sample sheet" not in content
	assert "Open ResPro Report" not in content

	(stage / "sample_plot_ref.svg").write_text("<svg></svg>")
	(stage / "sample_plot_ref_bed_data_1.tabular").write_text(BED)
	stage2, variants2 = _write_sample(tmp_path, "sample_2")
	render_report(
		{"sample": stage, "sample_2": stage2},
		{"sample": variants, "sample_2": variants2},
		{"sample": ["ref"], "sample_2": ["ref"]},
		True,
		output,
		respro_reports={"sample_2": "sample_2/sample_2.respro.html"},
	)
	content = output.read_text()
	assert "select-sample" in content
	assert 'data-sample="sample_2" hidden' in content
	assert 'href="sample_2/sample_2.respro.html" target="_blank" rel="noopener noreferrer" hidden>Open ResPro Report</a>' in content
	assert content.count("Open ResPro Report") == 1
	assert "data:image/svg+xml;base64," in content
	assert "amp_1" in content
	assert content.count("Amplicon recovery") == 1


def test_cli_main_renders_report(tmp_path: Path) -> None:
	output = tmp_path / "cli_report.html"
	stage, variants = _write_sample(tmp_path)
	(stage / "sample_plot_ref.svg").write_text("<svg></svg>")

	bamdash_txt = tmp_path / "bamdash_dirs.txt"
	variants_txt = tmp_path / "variants.txt"
	references_txt = tmp_path / "references.txt"
	config_yaml = tmp_path / "report_config.yaml"
	sample_sheet_csv = tmp_path / "sample_sheet.csv"
	respro_txt = tmp_path / "respro_reports.txt"
	bamdash_txt.write_text(f"{stage}\n")
	variants_txt.write_text(f"{variants}\n")
	references_txt.write_text(json.dumps({"sample": ["ref"]}))
	config_yaml.write_text("res_dir: results\nbamdash_coverage: 10\nextra: '<tag>'\n")
	sample_sheet_csv.write_text("sample,reference_fasta,r1,r2\nsample,ref.fa,r1.fq,r2.fq\n")
	respro_txt.write_text(json.dumps({"sample": "sample/sample.respro.html"}))

	main([
		"--bamdash-dirs", str(bamdash_txt),
		"--variants", str(variants_txt),
		"--references", str(references_txt),
		"--respro-reports", str(respro_txt),
		"--config", str(config_yaml),
		"--sample-sheet", str(sample_sheet_csv),
		"--no-svg",
		"--output", str(output),
	])
	content = output.read_text()
	assert "181" in content and "42" in content
	assert "data:image/svg+xml;base64," not in content
	assert "<details class=\"config-block\">" in content
	assert "res_dir: results\nbamdash_coverage: 10\nextra: &#39;&lt;tag&gt;&#39;" in content
	assert "Sample sheet" in content
	assert "sample,reference_fasta,r1,r2\nsample,ref.fa,r1.fq,r2.fq" in content
	assert 'href="sample/sample.respro.html" target="_blank" rel="noopener noreferrer">Open ResPro Report</a>' in content
	assert "Generated at: " in content
