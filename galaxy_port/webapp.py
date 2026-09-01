"""FastAPI web application for preparing Snakemake runs."""

from __future__ import annotations

import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

import yaml
from fastapi import FastAPI, File, HTTPException, Request, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse, RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel, Field

from galaxy_port.foundation import (
	AnalysisDefinition,
	PipelineConfigResolver,
	PrimerSelection,
	UploadedSample,
	build_run_config,
	create_run_directory,
	generate_snakemake_command,
	load_reference_catalog,
	load_tool_catalog,
	match_fastqs_to_samples,
	parse_illumina_samplesheet,
	write_run_config,
	write_samples_csv,
)


class AppSettings(BaseModel):
	"""Configuration for the run-preparation web application."""

	storage_root: Path = Path("/tmp/galaxy_port_runs")
	tus_upload_root: Path = Path("/tmp/galaxy_port_tusd")
	base_config_path: Path = Path("config.yaml")
	reference_catalog_path: Path = Path("galaxy_port/catalogs/references.yaml")
	tool_catalog_path: Path = Path("galaxy_port/catalogs/tool_options.yaml")
	repo_root: Path = Path(__file__).resolve().parents[1]
	tus_hook_secret: str | None = None


def load_settings_from_env() -> AppSettings:
	"""Build application settings from environment variables."""

	return AppSettings(
		storage_root=Path(os.environ.get("GALAXY_PORT_STORAGE_ROOT", "/tmp/galaxy_port_runs")),
		tus_upload_root=Path(os.environ.get("GALAXY_PORT_TUSD_UPLOAD_ROOT", "/tmp/galaxy_port_tusd")),
		base_config_path=Path(os.environ.get("GALAXY_PORT_BASE_CONFIG", "config.yaml")),
		reference_catalog_path=Path(os.environ.get("GALAXY_PORT_REFERENCE_CATALOG", "galaxy_port/catalogs/references.yaml")),
		tool_catalog_path=Path(os.environ.get("GALAXY_PORT_TOOL_CATALOG", "galaxy_port/catalogs/tool_options.yaml")),
		repo_root=Path(os.environ.get("GALAXY_PORT_REPO_ROOT", str(Path(__file__).resolve().parents[1]))),
		tus_hook_secret=os.environ.get("GALAXY_PORT_TUSD_HOOK_SECRET") or None,
	)


class PlannedFile(BaseModel):
	"""A selected FASTQ candidate from the browser before upload."""

	name: str
	size: int


class UploadPlanRequest(BaseModel):
	"""Request payload for FASTQ pre-flight validation."""

	files: list[PlannedFile]


class AnalysisRequest(BaseModel):
	"""One user-configured analysis row."""

	analysis_id: str
	uploaded_sample_id: str
	reference_id: str
	primer_mode: Literal["trim", "none"]
	primer_scheme_id: str | None = None
	tool_overrides: dict[str, dict[str, Any]] = Field(default_factory=dict)


class PrepareRunRequest(BaseModel):
	"""Payload used to finalize and prepare a run."""

	analyses: list[AnalysisRequest]
	run_overrides: dict[str, dict[str, Any]] = Field(default_factory=dict)


class TusUploadInfo(BaseModel):
	"""Normalized minimal tus upload information extracted from hook payloads."""

	upload_id: str
	run_id: str
	filename: str
	size: int
	storage_path: Path | None = None


def create_app(settings: AppSettings | None = None) -> FastAPI:
	"""Create and configure the FastAPI application instance."""

	effective_settings = settings or AppSettings()
	templates = Jinja2Templates(directory=str(Path(__file__).resolve().parent / "web" / "templates"))
	app = FastAPI(title="Galaxy Port Run Preparation")
	app.state.settings = effective_settings
	app.state.templates = templates

	static_dir = Path(__file__).resolve().parent / "web" / "static"
	app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")

	vendor_dir = effective_settings.repo_root / "node_modules"
	if vendor_dir.exists():
		app.mount("/vendor", StaticFiles(directory=str(vendor_dir)), name="vendor")

	@app.get("/", response_class=HTMLResponse)
	def index(request: Request) -> HTMLResponse:
		runs_root = effective_settings.storage_root / "runs"
		runs_root.mkdir(parents=True, exist_ok=True)
		runs = sorted([entry.name for entry in runs_root.glob("*") if entry.is_dir()], reverse=True)
		return templates.TemplateResponse(
			request,
			"index.html",
			{
				"request": request,
				"runs": runs,
			},
		)

	@app.post("/runs")
	def create_run() -> RedirectResponse:
		paths = create_run_directory(effective_settings.storage_root)
		state = _load_run_state(paths.run_json)
		state["sample_ids"] = []
		state["upload_plan"] = {}
		state["completed_uploads"] = {}
		state["warnings"] = []
		_save_run_state(paths.run_json, state)
		return RedirectResponse(url=f"/runs/{paths.run_id}", status_code=303)

	@app.get("/runs/{run_id}", response_class=HTMLResponse)
	def run_detail(run_id: str, request: Request) -> HTMLResponse:
		run_dir = _run_dir(effective_settings, run_id)
		state = _load_run_state(_run_json_path(run_dir))
		reference_catalog = load_reference_catalog(effective_settings.reference_catalog_path)
		tool_catalog = load_tool_catalog(effective_settings.tool_catalog_path)
		context = {
			"request": request,
			"run_id": run_id,
			"state": state,
			"references": _reference_options(reference_catalog),
			"tool_catalog": _tool_catalog_for_template(tool_catalog),
			"plan_json": json.dumps(state.get("upload_plan") or {}),
		}
		return templates.TemplateResponse(request, "run_detail.html", context)

	@app.post("/runs/{run_id}/samplesheet")
	async def upload_samplesheet(run_id: str, sample_sheet: UploadFile = File(...)) -> RedirectResponse:
		run_dir = _run_dir(effective_settings, run_id)
		run_json_path = _run_json_path(run_dir)
		state = _load_run_state(run_json_path)
		if state.get("status") != "staging":
			raise HTTPException(status_code=409, detail="Run is not in staging state")
		if not sample_sheet.filename:
			raise HTTPException(status_code=400, detail="Missing filename")
		if not sample_sheet.filename.lower().endswith(".csv"):
			raise HTTPException(status_code=400, detail="SampleSheet must be a CSV file")
		input_dir = run_dir / "input"
		input_dir.mkdir(parents=True, exist_ok=True)
		target = input_dir / "SampleSheet.csv"
		content = await sample_sheet.read()
		target.write_bytes(content)
		sample_ids = parse_illumina_samplesheet(target)
		state["sample_ids"] = sample_ids
		state["upload_plan"] = {}
		state["completed_uploads"] = {}
		state["warnings"] = []
		_save_run_state(run_json_path, state)
		return RedirectResponse(url=f"/runs/{run_id}", status_code=303)

	@app.post("/runs/{run_id}/upload-plan")
	def create_upload_plan(run_id: str, payload: UploadPlanRequest) -> JSONResponse:
		run_dir = _run_dir(effective_settings, run_id)
		run_json_path = _run_json_path(run_dir)
		state = _load_run_state(run_json_path)
		if state.get("status") != "staging":
			raise HTTPException(status_code=409, detail="Run is not in staging state")
		sample_ids = state.get("sample_ids") or []
		if not sample_ids:
			raise HTTPException(status_code=400, detail="Upload a valid SampleSheet first")
		plan = _validate_upload_plan(
			effective_settings,
			run_id,
			sample_ids,
			payload.files,
		)
		uploads_dir = run_dir / ".uploads"
		uploads_dir.mkdir(parents=True, exist_ok=True)
		(uploads_dir / "upload-plan.json").write_text(json.dumps(plan, indent=2))
		state["upload_plan"] = plan
		state["warnings"] = plan.get("warnings", [])
		_save_run_state(run_json_path, state)
		return JSONResponse(plan)

	@app.get("/runs/{run_id}/upload-status-fragment", response_class=HTMLResponse)
	def upload_status_fragment(run_id: str, request: Request) -> HTMLResponse:
		run_dir = _run_dir(effective_settings, run_id)
		state = _load_run_state(_run_json_path(run_dir))
		completed_uploads = state.get("completed_uploads") or {}
		plan = state.get("upload_plan") or {}
		required_files = set((plan.get("required_files") or {}).keys())
		confirmed = sorted(name for name in completed_uploads if name in required_files)
		missing = sorted(required_files - set(confirmed))
		return templates.TemplateResponse(
			request,
			"_upload_status.html",
			{
				"request": request,
				"confirmed": confirmed,
				"missing": missing,
			},
		)

	@app.get("/runs/{run_id}/analysis-row", response_class=HTMLResponse)
	def analysis_row_fragment(run_id: str, request: Request, index: int = 0) -> HTMLResponse:
		run_dir = _run_dir(effective_settings, run_id)
		state = _load_run_state(_run_json_path(run_dir))
		reference_catalog = load_reference_catalog(effective_settings.reference_catalog_path)
		tool_catalog = load_tool_catalog(effective_settings.tool_catalog_path)
		return templates.TemplateResponse(
			request,
			"_analysis_row.html",
			{
				"request": request,
				"index": index,
				"sample_ids": state.get("sample_ids") or [],
				"references": _reference_options(reference_catalog),
				"tool_catalog": _tool_catalog_for_template(tool_catalog),
			},
		)

	@app.post("/internal/tusd/hooks/pre-create")
	async def tusd_pre_create(request: Request) -> JSONResponse:
		_hook_guard(request, effective_settings)
		payload = await request.json()
		info = _extract_tusd_upload_info(payload)
		state, run_dir = _load_run_state_for_hook(effective_settings, info.run_id)
		plan = state.get("upload_plan") or {}
		required_files = plan.get("required_files") or {}
		completed_uploads = state.get("completed_uploads") or {}
		if state.get("status") != "staging":
			return _reject_upload("Run is not in staging state")
		if not _is_safe_filename(info.filename):
			return _reject_upload("Unsafe upload filename")
		if info.filename not in required_files:
			return _reject_upload("Filename is not part of approved upload plan")
		expected_size = int(required_files[info.filename])
		if expected_size != info.size:
			return _reject_upload("Upload size does not match approved upload plan")
		if info.filename in completed_uploads:
			return _reject_upload("File already completed")
		if (run_dir / "input" / info.filename).exists():
			return _reject_upload("Destination filename already exists")
		return JSONResponse({"ok": True})

	@app.post("/internal/tusd/hooks/post-finish")
	async def tusd_post_finish(request: Request) -> JSONResponse:
		_hook_guard(request, effective_settings)
		payload = await request.json()
		info = _extract_tusd_upload_info(payload)
		state, run_dir = _load_run_state_for_hook(effective_settings, info.run_id)
		plan = state.get("upload_plan") or {}
		required_files = plan.get("required_files") or {}
		if info.filename not in required_files:
			raise HTTPException(status_code=400, detail="Filename is not part of approved upload plan")
		if int(required_files[info.filename]) != info.size:
			raise HTTPException(status_code=400, detail="Uploaded size does not match plan")
		target = run_dir / "input" / info.filename
		target.parent.mkdir(parents=True, exist_ok=True)
		if info.storage_path and info.storage_path.exists() and not target.exists():
			info.storage_path.replace(target)
		if target.exists() and target.stat().st_size != info.size:
			raise HTTPException(status_code=400, detail="Stored file size does not match hook size")
		if not target.exists():
			raise HTTPException(status_code=400, detail="Upload hook payload did not include a usable storage path")
		completed_uploads = state.get("completed_uploads") or {}
		completed_uploads[info.filename] = {
			"upload_id": info.upload_id,
			"filename": info.filename,
			"size": info.size,
			"path": str(target),
			"completed_at": datetime.now(timezone.utc).isoformat(),
		}
		state["completed_uploads"] = completed_uploads
		_save_run_state(_run_json_path(run_dir), state)
		upload_record_path = run_dir / ".uploads" / f"{info.upload_id}.json"
		upload_record_path.write_text(json.dumps(completed_uploads[info.filename], indent=2))
		return JSONResponse({"ok": True})

	@app.post("/runs/{run_id}/prepare")
	def prepare_run(run_id: str, payload: PrepareRunRequest) -> JSONResponse:
		run_dir = _run_dir(effective_settings, run_id)
		run_json_path = _run_json_path(run_dir)
		state = _load_run_state(run_json_path)
		if state.get("status") != "staging":
			raise HTTPException(status_code=409, detail="Run is not in staging state")
		plan = state.get("upload_plan") or {}
		required_files = plan.get("required_files") or {}
		completed_uploads = state.get("completed_uploads") or {}
		missing_files = sorted(set(required_files.keys()) - set(completed_uploads.keys()))
		if missing_files:
			raise HTTPException(status_code=400, detail=f"Run cannot be prepared. Missing completed uploads: {missing_files}")
		uploaded_samples = _uploaded_samples_from_plan(plan, run_dir)
		if not payload.analyses:
			raise HTTPException(status_code=400, detail="At least one analysis is required")
		analysis_models: list[AnalysisDefinition] = []
		for analysis in payload.analyses:
			primer = PrimerSelection(mode=analysis.primer_mode, primer_scheme_id=analysis.primer_scheme_id)
			analysis_models.append(
				AnalysisDefinition(
					analysis_id=analysis.analysis_id,
					uploaded_sample_id=analysis.uploaded_sample_id,
					reference_id=analysis.reference_id,
					primer=primer,
					tool_overrides=analysis.tool_overrides,
				)
			)
		reference_catalog = load_reference_catalog(effective_settings.reference_catalog_path)
		write_samples_csv(run_dir / "samples.csv", analysis_models, uploaded_samples, reference_catalog)
		base_config = yaml.safe_load(effective_settings.base_config_path.read_text()) or {}
		prepared_config = build_run_config(
			base_config=base_config,
			sample_sheet_path=run_dir / "samples.csv",
			res_dir=run_dir / "results",
			run_overrides=payload.run_overrides,
			analyses=analysis_models,
			tool_catalog=load_tool_catalog(effective_settings.tool_catalog_path),
		)
		write_run_config(run_dir / "config.yaml", prepared_config)
		PipelineConfigResolver(prepared_config, tool_catalog=load_tool_catalog(effective_settings.tool_catalog_path))
		command = generate_snakemake_command(
			repo_root=effective_settings.repo_root,
			config_path=run_dir / "config.yaml",
			singularity_bind_paths=[run_dir, effective_settings.repo_root],
		)
		state["status"] = "prepared"
		state["analyses"] = [analysis.model_dump() for analysis in analysis_models]
		state["run_overrides"] = payload.run_overrides
		state["snakemake_command"] = command
		_save_run_state(run_json_path, state)
		return JSONResponse(
			{
				"status": "prepared",
				"run_id": run_id,
				"snakemake_command": command,
				"samples_csv": str(run_dir / "samples.csv"),
				"config_yaml": str(run_dir / "config.yaml"),
			}
		)

	return app


def _uploaded_samples_from_plan(plan: dict[str, Any], run_dir: Path) -> list[UploadedSample]:
	uploaded_samples: list[UploadedSample] = []
	for entry in plan.get("uploaded_samples") or []:
		uploaded_samples.append(
			UploadedSample(
				sample_id=entry["sample_id"],
				r1=run_dir / "input" / entry["r1"],
				r2=run_dir / "input" / entry["r2"],
			)
		)
	return uploaded_samples


def _validate_upload_plan(
	settings: AppSettings,
	run_id: str,
	sample_ids: list[str],
	selected_files: list[PlannedFile],
) -> dict[str, Any]:
	if not selected_files:
		raise HTTPException(status_code=400, detail="Select at least one FASTQ file")
	seen_names: set[str] = set()
	for file_info in selected_files:
		if file_info.name in seen_names:
			raise HTTPException(status_code=400, detail=f"Duplicate filename selected: {file_info.name}")
		seen_names.add(file_info.name)
		if file_info.size <= 0:
			raise HTTPException(status_code=400, detail=f"Zero-byte file selected: {file_info.name}")
		if not file_info.name.lower().endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq")):
			raise HTTPException(status_code=400, detail=f"Unsupported FASTQ suffix: {file_info.name}")
	try:
		match_result = match_fastqs_to_samples(sample_ids, [Path(item.name) for item in selected_files])
	except ValueError as exc:
		raise HTTPException(status_code=400, detail=str(exc)) from exc
	file_size_map = {item.name: int(item.size) for item in selected_files}
	required_files: dict[str, int] = {}
	uploaded_samples: list[dict[str, Any]] = []
	for sample in match_result.uploaded_samples:
		r1_name = Path(sample.r1).name
		r2_name = Path(sample.r2).name
		required_files[r1_name] = file_size_map[r1_name]
		required_files[r2_name] = file_size_map[r2_name]
		uploaded_samples.append(
			{
				"sample_id": sample.sample_id,
				"r1": r1_name,
				"r2": r2_name,
			}
		)
	total_expected_size = sum(required_files.values())
	free_space = shutil.disk_usage(settings.storage_root).free
	if free_space < total_expected_size:
		raise HTTPException(status_code=400, detail="Not enough free storage space for required FASTQ uploads")
	return {
		"run_id": run_id,
		"required_files": required_files,
		"uploaded_samples": uploaded_samples,
		"warnings": [f"Extra FASTQs ignored: {match_result.extra_fastqs}"] if match_result.extra_fastqs else [],
		"total_expected_size": total_expected_size,
		"free_space": free_space,
	}


def _run_dir(settings: AppSettings, run_id: str) -> Path:
	run_dir = settings.storage_root / "runs" / run_id
	if not run_dir.exists():
		raise HTTPException(status_code=404, detail=f"Unknown run_id '{run_id}'")
	return run_dir


def _run_json_path(run_dir: Path) -> Path:
	return run_dir / "run.json"


def _load_run_state(run_json_path: Path) -> dict[str, Any]:
	if not run_json_path.exists():
		raise HTTPException(status_code=404, detail="Run metadata not found")
	return json.loads(run_json_path.read_text())


def _save_run_state(run_json_path: Path, state: dict[str, Any]) -> None:
	run_json_path.write_text(json.dumps(state, indent=2, sort_keys=True))


def _reference_options(reference_catalog: Any) -> list[dict[str, Any]]:
	options: list[dict[str, Any]] = []
	for reference_id, reference in reference_catalog.references.items():
		options.append(
			{
				"reference_id": reference_id,
				"name": reference.name,
				"primer_schemes": [
					{
						"primer_scheme_id": primer_scheme_id,
						"name": primer.name,
					}
					for primer_scheme_id, primer in reference.primer_schemes.items()
				],
			}
		)
	return options


def _tool_catalog_for_template(tool_catalog: Any) -> list[dict[str, Any]]:
	tools: list[dict[str, Any]] = []
	for tool_id, tool in tool_catalog.tools.items():
		tools.append(
			{
				"tool_id": tool_id,
				"name": tool.name,
				"category": tool.category,
				"options": [
					{
						"option_id": option_id,
						"label": option.label,
						"type": option.type,
						"default": option.default,
						"minimum": option.minimum,
						"maximum": option.maximum,
						"choices": option.choices,
					}
					for option_id, option in tool.options.items()
				],
			}
		)
	return tools


def _hook_guard(request: Request, settings: AppSettings) -> None:
	if settings.tus_hook_secret is None:
		return
	secret = request.headers.get("x-galaxy-port-hook-secret", "")
	if secret != settings.tus_hook_secret:
		raise HTTPException(status_code=403, detail="Invalid hook secret")


def _load_run_state_for_hook(settings: AppSettings, run_id: str) -> tuple[dict[str, Any], Path]:
	run_dir = settings.storage_root / "runs" / run_id
	run_json_path = run_dir / "run.json"
	if not run_json_path.exists():
		raise HTTPException(status_code=404, detail="Unknown run ID")
	return _load_run_state(run_json_path), run_dir


def _extract_tusd_upload_info(payload: dict[str, Any]) -> TusUploadInfo:
	upload = payload.get("Upload") or payload.get("upload") or {}
	metadata = upload.get("MetaData") or upload.get("metadata") or upload.get("metaData") or {}
	filename = metadata.get("filename", "")
	run_id = metadata.get("run_id", "")
	upload_id = str(upload.get("ID") or upload.get("id") or "")
	size = int(upload.get("Size") or upload.get("size") or 0)
	storage = upload.get("Storage") or upload.get("storage") or {}
	storage_path_text = storage.get("Path") or storage.get("path")
	if not upload_id or not filename or not run_id:
		raise HTTPException(status_code=400, detail="Missing upload metadata for hook event")
	storage_path = Path(storage_path_text) if storage_path_text else None
	return TusUploadInfo(
		upload_id=upload_id,
		run_id=run_id,
		filename=filename,
		size=size,
		storage_path=storage_path,
	)


def _reject_upload(message: str) -> JSONResponse:
	return JSONResponse(
		{
			"RejectUpload": True,
			"HTTPResponse": {
				"StatusCode": 400,
				"Body": message,
			},
		},
		status_code=200,
	)


def _is_safe_filename(name: str) -> bool:
	if not name or name.strip() != name:
		return False
	if "/" in name or "\\" in name:
		return False
	if ".." in name:
		return False
	path = Path(name)
	if path.is_absolute():
		return False
	return True


app = create_app(load_settings_from_env())