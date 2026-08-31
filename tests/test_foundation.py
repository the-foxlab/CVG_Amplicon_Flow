from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import pytest
import yaml

from galaxy_port.foundation import (
	AnalysisDefinition,
	PipelineConfigResolver,
	PrimerSelection,
	RunDefinition,
	build_pipeline_samples_map,
	build_run_config,
	create_run_directory,
	generate_samples_dataframe,
	generate_snakemake_command,
	load_reference_catalog,
	match_fastqs_to_samples,
	parse_illumina_samplesheet,
	write_run_config,
	write_samples_csv,
	UploadedSample,
)


def test_parse_illumina_samplesheet_reads_data_section(tmp_path: Path) -> None:
	sample_sheet = tmp_path / "SampleSheet.csv"
	sample_sheet.write_text(
		"[Header]\n"
		"IEMFileVersion,4\n"
		"[Data]\n"
		"Sample_ID,Sample_Name\n"
		"alpha,Alpha\n"
		"beta,Beta\n"
	)

	assert parse_illumina_samplesheet(sample_sheet) == ["alpha", "beta"]


def test_parse_illumina_samplesheet_rejects_duplicate_ids(tmp_path: Path) -> None:
	sample_sheet = tmp_path / "SampleSheet.csv"
	sample_sheet.write_text("Sample_ID\nalpha\nalpha\n")

	with pytest.raises(ValueError, match="duplicate sample IDs"):
		parse_illumina_samplesheet(sample_sheet)


def test_match_fastqs_to_samples_reports_extra_fastq(tmp_path: Path) -> None:
	paths = [
		tmp_path / "alpha_S1_L001_R1_001.fastq.gz",
		tmp_path / "alpha_S1_L001_R2_001.fastq.gz",
		tmp_path / "extra_S1_L001_R1_001.fastq.gz",
	]
	result = match_fastqs_to_samples(["alpha"], paths)

	assert [sample.sample_id for sample in result.uploaded_samples] == ["alpha"]
	assert result.extra_fastqs == ["extra_S1_L001_R1_001.fastq.gz"]


def test_match_fastqs_to_samples_rejects_missing_pair(tmp_path: Path) -> None:
	with pytest.raises(ValueError, match="Missing R2 FASTQ"):
		match_fastqs_to_samples(["alpha"], [tmp_path / "alpha_S1_L001_R1_001.fastq.gz"])


def test_run_definition_rejects_duplicate_analysis_ids(tmp_path: Path) -> None:
	uploaded_sample = UploadedSample(
		sample_id="alpha",
		r1=tmp_path / "alpha_S1_L001_R1_001.fastq.gz",
		r2=tmp_path / "alpha_S1_L001_R2_001.fastq.gz",
	)
	analysis = AnalysisDefinition(
		analysis_id="alpha_trimmed",
		uploaded_sample_id="alpha",
		reference_id="hsv2_both",
		primer=PrimerSelection(mode="none"),
	)

	with pytest.raises(ValueError, match="unique analysis_id"):
		RunDefinition(uploaded_samples=[uploaded_sample], analyses=[analysis, analysis])


def test_pipeline_config_resolver_uses_override_precedence_and_legacy_fallback() -> None:
	config = {
		"min_base_quality": 20,
		"bwa_mem_extra": "-Y",
		"run_overrides": {
			"fastp": {"qualified_quality_phred": 30},
			"bwa_mem": {"band_width": 80},
		},
		"analysis_overrides": {
			"alpha_trimmed": {
				"fastp": {"qualified_quality_phred": 40},
				"bwa_mem": {"min_seed_length": 17},
			}
		},
	}
	resolver = PipelineConfigResolver(config)

	assert resolver.get_analysis_parameter("alpha_trimmed", "fastp", "qualified_quality_phred") == 40
	assert resolver.get_analysis_parameter("alpha_other", "fastp", "qualified_quality_phred") == 30
	assert resolver.get_analysis_parameter("alpha_other", "ivar_consensus", "min_depth") == 20
	assert resolver.render_tool_args("alpha_trimmed", "bwa_mem", "bwa_mem_extra").endswith("-Y")
	assert "-k 17" in resolver.render_tool_args("alpha_trimmed", "bwa_mem", "bwa_mem_extra")
	assert "-w 80" in resolver.render_tool_args("alpha_trimmed", "bwa_mem", "bwa_mem_extra")


def test_pipeline_config_resolver_rejects_unknown_options() -> None:
	with pytest.raises(ValueError, match="Unknown option"):
		PipelineConfigResolver(
			{
				"run_overrides": {
					"fastp": {"raw_args": "-x"},
				}
			}
		)


def test_build_pipeline_samples_map_validates_fasta_and_bed(tmp_path: Path) -> None:
	fasta = tmp_path / "ref.fasta"
	bed = tmp_path / "primers.bed"
	sample_csv = tmp_path / "samples.csv"
	fasta.write_text(">refA\nACGT\n>refB\nTGCA\n")
	bed.write_text("refA\t0\t4\n")
	sample_csv.write_text(
		"sample,reference_fasta,bed,r1,r2\n"
		f"alpha,{fasta},{bed},{tmp_path / 'alpha_R1.fastq.gz'},{tmp_path / 'alpha_R2.fastq.gz'}\n"
	)

	samples_map = build_pipeline_samples_map(sample_csv)

	assert list(samples_map) == ["alpha"]
	assert samples_map["alpha"]["references"] == ["refA"]
	assert samples_map["alpha"]["bed"] == str(bed)


def test_generate_artifacts_and_command(tmp_path: Path) -> None:
	reference_catalog_path = tmp_path / "references.yaml"
	reference_catalog_path.write_text(
		yaml.safe_dump(
			{
				"references": {
					"hsv2_both": {
						"name": "HSV-2 both",
						"fasta": str(tmp_path / "HSV_2_both.fasta"),
						"primer_schemes": {
							"hsv2_both": {
								"name": "HSV-2 both primers",
								"bed": str(tmp_path / "HSV_2_both.bed"),
							}
						},
					}
				}
			},
			sort_keys=False,
		)
	)
	(tmp_path / "HSV_2_both.fasta").write_text(">refA\nACGT\n")
	(tmp_path / "HSV_2_both.bed").write_text("refA\t0\t4\n")
	reference_catalog = load_reference_catalog(reference_catalog_path)
	uploaded_sample = UploadedSample(
		sample_id="1545554-HSV2_S82",
		r1=tmp_path / "1545554-HSV2_S82_L001_R1_001.fastq.gz",
		r2=tmp_path / "1545554-HSV2_S82_L001_R2_001.fastq.gz",
	)
	analyses = [
		AnalysisDefinition(
			analysis_id="1545554-HSV2_S82_trimmed",
			uploaded_sample_id="1545554-HSV2_S82",
			reference_id="hsv2_both",
			primer=PrimerSelection(mode="trim", primer_scheme_id="hsv2_both"),
			tool_overrides={"ivar_consensus": {"min_depth": 50}},
		),
		AnalysisDefinition(
			analysis_id="1545554-HSV2_S82_not_trimmed",
			uploaded_sample_id="1545554-HSV2_S82",
			reference_id="hsv2_both",
			primer=PrimerSelection(mode="none"),
		),
	]
	dataframe = generate_samples_dataframe(analyses, [uploaded_sample], reference_catalog)
	assert dataframe["sample"].tolist() == [
		"1545554-HSV2_S82_trimmed",
		"1545554-HSV2_S82_not_trimmed",
	]
	assert dataframe["bed"].tolist()[1] == ""
	run_paths = create_run_directory(
		tmp_path,
		run_id="20260831-122541-a1b2c3",
		created_at=datetime(2026, 8, 31, 12, 25, 41, tzinfo=timezone.utc),
	)
	write_samples_csv(run_paths.samples_csv, analyses, [uploaded_sample], reference_catalog)
	config = build_run_config(
		base_config={"min_base_quality": 20, "min_depth": 20},
		sample_sheet_path=run_paths.samples_csv,
		res_dir=run_paths.results_dir,
		run_overrides={"fastp": {"qualified_quality_phred": 25}},
		analyses=analyses,
	)
	write_run_config(run_paths.config_yaml, config)
	command = generate_snakemake_command(
		repo_root=Path.cwd(),
		config_path=run_paths.config_yaml,
		singularity_bind_paths=[tmp_path],
	)

	assert run_paths.run_json.exists()
	assert run_paths.samples_csv.exists()
	assert run_paths.config_yaml.exists()
	loaded_config = yaml.safe_load(run_paths.config_yaml.read_text())
	assert loaded_config["analysis_overrides"]["1545554-HSV2_S82_trimmed"]["ivar_consensus"]["min_depth"] == 50
	assert "--configfile" in command
	assert str(run_paths.config_yaml) in command