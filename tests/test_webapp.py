from __future__ import annotations

import json
from pathlib import Path

import yaml
from fastapi.testclient import TestClient

from galaxy_port.webapp import AppSettings, create_app


def _make_settings(tmp_path: Path) -> AppSettings:
	repo_root = Path(__file__).resolve().parents[1]
	storage_root = tmp_path / "runtime"
	tus_upload_root = tmp_path / "tusd"
	base_config_path = tmp_path / "base_config.yaml"
	reference_catalog_path = tmp_path / "references.yaml"
	tool_catalog_path = repo_root / "galaxy_port" / "catalogs" / "tool_options.yaml"
	(storage_root / "runs").mkdir(parents=True, exist_ok=True)
	tus_upload_root.mkdir(parents=True, exist_ok=True)
	base_config_path.write_text(
		yaml.safe_dump(
			{
				"sample_sheet": "sample_sheet.csv",
				"res_dir": "/tmp/results",
			},
			sort_keys=False,
		)
	)
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
	return AppSettings(
		storage_root=storage_root,
		tus_upload_root=tus_upload_root,
		base_config_path=base_config_path,
		reference_catalog_path=reference_catalog_path,
		tool_catalog_path=tool_catalog_path,
		repo_root=repo_root,
		tus_hook_secret="test-secret",
	)


def _create_run(client: TestClient) -> str:
	response = client.post("/runs", follow_redirects=False)
	assert response.status_code == 303
	location = response.headers["location"]
	return location.rsplit("/", 1)[-1]


def _upload_samplesheet(client: TestClient, run_id: str) -> None:
	sample_sheet = (
		"[Header]\n"
		"IEMFileVersion,4\n"
		"[Data]\n"
		"Sample_ID,Sample_Name\n"
		"1545554-HSV2_S82,SampleA\n"
	)
	response = client.post(
		f"/runs/{run_id}/samplesheet",
		files={"sample_sheet": ("SampleSheet.csv", sample_sheet, "text/csv")},
		follow_redirects=False,
	)
	assert response.status_code == 303


def _create_plan(client: TestClient, run_id: str) -> dict:
	payload = {
		"files": [
			{"name": "1545554-HSV2_S82_S1_L001_R1_001.fastq.gz", "size": 4},
			{"name": "1545554-HSV2_S82_S1_L001_R2_001.fastq.gz", "size": 4},
		]
	}
	response = client.post(f"/runs/{run_id}/upload-plan", json=payload)
	assert response.status_code == 200
	plan = response.json()
	assert sorted(plan["required_files"].keys()) == sorted(
		[
			"1545554-HSV2_S82_S1_L001_R1_001.fastq.gz",
			"1545554-HSV2_S82_S1_L001_R2_001.fastq.gz",
		]
	)
	return plan


def _post_hooks(client: TestClient, settings: AppSettings, run_id: str, filename: str, upload_id: str) -> None:
	upload_path = settings.tus_upload_root / upload_id
	upload_path.write_bytes(b"ACGT")
	headers = {"x-galaxy-port-hook-secret": "test-secret"}
	pre_create = {
		"Upload": {
			"ID": upload_id,
			"Size": 4,
			"MetaData": {
				"filename": filename,
				"run_id": run_id,
			},
		}
	}
	response_pre = client.post("/internal/tusd/hooks/pre-create", headers=headers, json=pre_create)
	assert response_pre.status_code == 200
	assert response_pre.json().get("ok") is True
	post_finish = {
		"Upload": {
			"ID": upload_id,
			"Size": 4,
			"MetaData": {
				"filename": filename,
				"run_id": run_id,
			},
			"Storage": {
				"Path": str(upload_path),
			},
		}
	}
	response_post = client.post("/internal/tusd/hooks/post-finish", headers=headers, json=post_finish)
	assert response_post.status_code == 200
	assert response_post.json().get("ok") is True


def test_prepare_run_full_flow(tmp_path: Path) -> None:
	settings = _make_settings(tmp_path)
	client = TestClient(create_app(settings))
	run_id = _create_run(client)
	_upload_samplesheet(client, run_id)
	plan = _create_plan(client, run_id)
	assert len(plan["uploaded_samples"]) == 1
	_post_hooks(client, settings, run_id, "1545554-HSV2_S82_S1_L001_R1_001.fastq.gz", "upload-r1")
	_post_hooks(client, settings, run_id, "1545554-HSV2_S82_S1_L001_R2_001.fastq.gz", "upload-r2")
	prepare_payload = {
		"analyses": [
			{
				"analysis_id": "1545554-HSV2_S82_trimmed",
				"uploaded_sample_id": "1545554-HSV2_S82",
				"reference_id": "hsv2_both",
				"primer_mode": "trim",
				"primer_scheme_id": "hsv2_both",
				"tool_overrides": {
					"ivar_consensus": {
						"min_depth": 50,
					}
				},
			},
			{
				"analysis_id": "1545554-HSV2_S82_not_trimmed",
				"uploaded_sample_id": "1545554-HSV2_S82",
				"reference_id": "hsv2_both",
				"primer_mode": "none",
				"primer_scheme_id": None,
				"tool_overrides": {},
			},
		],
		"run_overrides": {
			"fastp": {
				"qualified_quality_phred": 25,
			}
		},
	}
	response = client.post(f"/runs/{run_id}/prepare", json=prepare_payload)
	assert response.status_code == 200
	data = response.json()
	assert data["status"] == "prepared"
	run_dir = settings.storage_root / "runs" / run_id
	assert (run_dir / "samples.csv").exists()
	assert (run_dir / "config.yaml").exists()
	state = json.loads((run_dir / "run.json").read_text())
	assert state["status"] == "prepared"
	assert "snakemake" in state["snakemake_command"]


def test_tusd_precreate_rejects_unplanned_filename(tmp_path: Path) -> None:
	settings = _make_settings(tmp_path)
	client = TestClient(create_app(settings))
	run_id = _create_run(client)
	_upload_samplesheet(client, run_id)
	_create_plan(client, run_id)
	headers = {"x-galaxy-port-hook-secret": "test-secret"}
	payload = {
		"Upload": {
			"ID": "bad-1",
			"Size": 4,
			"MetaData": {
				"filename": "not_planned_S1_L001_R1_001.fastq.gz",
				"run_id": run_id,
			},
		}
	}
	response = client.post("/internal/tusd/hooks/pre-create", headers=headers, json=payload)
	assert response.status_code == 200
	assert response.json()["RejectUpload"] is True


def test_analysis_row_and_upload_status_fragments_render(tmp_path: Path) -> None:
	settings = _make_settings(tmp_path)
	client = TestClient(create_app(settings))
	run_id = _create_run(client)
	_upload_samplesheet(client, run_id)
	response_row = client.get(f"/runs/{run_id}/analysis-row?index=0")
	assert response_row.status_code == 200
	assert "Analysis 1" in response_row.text
	response_status = client.get(f"/runs/{run_id}/upload-status-fragment")
	assert response_status.status_code == 200
	assert "Upload Completion" in response_status.text
