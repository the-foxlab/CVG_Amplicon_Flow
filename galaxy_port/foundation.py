"""Domain and configuration helpers for the GUI MVP foundation."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import shlex
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from secrets import token_hex
from typing import Any, Iterable, Literal, Mapping, Sequence

import pandas as pd
import yaml
from Bio import SeqIO
from pydantic import BaseModel, Field, field_validator, model_validator


CATALOGS_DIR = Path(__file__).resolve().parent / "catalogs"
DEFAULT_REFERENCE_CATALOG_PATH = CATALOGS_DIR / "references.yaml"
DEFAULT_TOOL_CATALOG_PATH = CATALOGS_DIR / "tool_options.yaml"
FASTQ_NAME_PATTERN = re.compile(
	r"^(?P<prefix>.+)_S\d+(?:_L\d{3})?_R(?P<read>[12])_001(?P<suffix>\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$"
)
ANALYSIS_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
MISSING_TEXT = {"", "nan", "none", "null"}


def _model_validate(model_class: type[BaseModel], data: Any) -> BaseModel:
	if hasattr(model_class, "model_validate"):
		return model_class.model_validate(data)
	return model_class.parse_obj(data)


def _model_dump(model: BaseModel) -> dict[str, Any]:
	if hasattr(model, "model_dump"):
		return model.model_dump()
	return model.dict()


def _is_missing(value: Any) -> bool:
	text = str(value).strip().lower()
	return text in MISSING_TEXT


def _normalize_path(value: Path | str) -> Path:
	return value if isinstance(value, Path) else Path(value)


def _safe_analysis_id(value: str) -> str:
	cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
	cleaned = cleaned.strip("._-")
	if not cleaned:
		raise ValueError("analysis_id cannot be empty after sanitization")
	if not ANALYSIS_ID_PATTERN.fullmatch(cleaned):
		raise ValueError(f"analysis_id '{value}' is not filesystem safe")
	return cleaned


class PrimerSelection(BaseModel):
	"""Primer-trimming choice for an analysis."""

	mode: Literal["trim", "none"]
	primer_scheme_id: str | None = None

	@model_validator(mode="after")
	def _validate_mode(self) -> "PrimerSelection":
		if self.mode == "trim" and not self.primer_scheme_id:
			raise ValueError("primer_scheme_id is required when mode='trim'")
		if self.mode == "none" and self.primer_scheme_id is not None:
			raise ValueError("primer_scheme_id must be omitted when mode='none'")
		return self


class UploadedSample(BaseModel):
	"""Uploaded FASTQ pair for one Illumina sample."""

	sample_id: str
	r1: Path
	r2: Path

	@field_validator("sample_id")
	@classmethod
	def _validate_sample_id(cls, value: str) -> str:
		value = value.strip()
		if not value:
			raise ValueError("sample_id cannot be empty")
		return value


class AnalysisDefinition(BaseModel):
	"""A single Snakemake analysis entry derived from an uploaded sample."""

	analysis_id: str
	uploaded_sample_id: str
	reference_id: str
	primer: PrimerSelection
	tool_overrides: dict[str, dict[str, Any]] = Field(default_factory=dict)

	@field_validator("analysis_id")
	@classmethod
	def _validate_analysis_id(cls, value: str) -> str:
		return _safe_analysis_id(value)


class PrimerSchemeDefinition(BaseModel):
	"""Primer scheme metadata from the repository catalog."""

	name: str
	bed: Path


class ReferenceDefinition(BaseModel):
	"""Reference metadata from the repository catalog."""

	name: str
	fasta: Path
	primer_schemes: dict[str, PrimerSchemeDefinition] = Field(default_factory=dict)


class ReferenceCatalog(BaseModel):
	"""Repository-controlled reference and primer catalog."""

	references: dict[str, ReferenceDefinition]

	def resolve_paths(self, reference_id: str, primer: PrimerSelection) -> tuple[Path, str]:
		"""Resolve the FASTA and BED path for one analysis."""

		if reference_id not in self.references:
			raise ValueError(f"Unknown reference_id '{reference_id}'")
		reference = self.references[reference_id]
		if primer.mode == "none":
			return reference.fasta, ""
		if primer.primer_scheme_id not in reference.primer_schemes:
			raise ValueError(
				f"Primer scheme '{primer.primer_scheme_id}' is not valid for reference '{reference_id}'"
			)
		return reference.fasta, str(reference.primer_schemes[primer.primer_scheme_id].bed)


class ToolOptionDefinition(BaseModel):
	"""Definition of one allowed tool option."""

	flag: str
	label: str
	type: Literal["bool", "int", "float", "string", "enum", "list"]
	default: Any = None
	minimum: float | None = None
	maximum: float | None = None
	choices: list[str] = Field(default_factory=list)
	help_text: str | None = None
	legacy_keys: list[str] = Field(default_factory=list)
	repeatable: bool = False
	item_type: Literal["bool", "int", "float", "string", "enum"] | None = None

	@model_validator(mode="after")
	def _validate_definition(self) -> "ToolOptionDefinition":
		option_type = self.type
		choices = self.choices or []
		item_type = self.item_type
		if option_type == "enum" and not choices:
			raise ValueError("enum options require at least one choice")
		if option_type != "enum" and choices:
			raise ValueError("choices are only valid for enum options")
		if option_type == "list" and item_type is None:
			raise ValueError("list options require item_type")
		if option_type != "list" and item_type is not None:
			raise ValueError("item_type is only valid for list options")
		return self

	def validate_value(self, value: Any, option_id: str) -> Any:
		"""Validate and normalize one option value."""

		if self.type == "bool":
			if not isinstance(value, bool):
				raise ValueError(f"Option '{option_id}' must be a boolean")
			return value
		if self.type == "int":
			if isinstance(value, bool) or not isinstance(value, int):
				raise ValueError(f"Option '{option_id}' must be an integer")
			self._validate_range(float(value), option_id)
			return value
		if self.type == "float":
			if isinstance(value, bool) or not isinstance(value, (int, float)):
				raise ValueError(f"Option '{option_id}' must be a float")
			normalized = float(value)
			self._validate_range(normalized, option_id)
			return normalized
		if self.type == "string":
			if not isinstance(value, str) or not value:
				raise ValueError(f"Option '{option_id}' must be a non-empty string")
			return value
		if self.type == "enum":
			if value not in self.choices:
				raise ValueError(
					f"Option '{option_id}' must be one of {self.choices}, got '{value}'"
				)
			return value
		if self.type == "list":
			if not isinstance(value, list):
				raise ValueError(f"Option '{option_id}' must be a list")
			return [self._validate_list_item(item, option_id) for item in value]
		raise ValueError(f"Unsupported option type '{self.type}'")

	def render_tokens(self, value: Any) -> list[str]:
		"""Render validated option values into CLI tokens."""

		if self.type == "bool":
			return [self.flag] if value else []
		if self.type == "list":
			tokens: list[str] = []
			for item in value:
				tokens.extend([self.flag, shlex.quote(str(item))])
			return tokens
		return [self.flag, shlex.quote(str(value))]

	def _validate_list_item(self, value: Any, option_id: str) -> Any:
		if self.item_type == "bool":
			if not isinstance(value, bool):
				raise ValueError(f"List option '{option_id}' items must be booleans")
			return value
		if self.item_type == "int":
			if isinstance(value, bool) or not isinstance(value, int):
				raise ValueError(f"List option '{option_id}' items must be integers")
			return value
		if self.item_type == "float":
			if isinstance(value, bool) or not isinstance(value, (int, float)):
				raise ValueError(f"List option '{option_id}' items must be floats")
			return float(value)
		if self.item_type == "string":
			if not isinstance(value, str) or not value:
				raise ValueError(f"List option '{option_id}' items must be non-empty strings")
			return value
		if self.item_type == "enum":
			if value not in self.choices:
				raise ValueError(
					f"List option '{option_id}' items must be one of {self.choices}, got '{value}'"
				)
			return value
		raise ValueError(f"Unsupported item_type '{self.item_type}'")

	def _validate_range(self, value: float, option_id: str) -> None:
		if self.minimum is not None and value < self.minimum:
			raise ValueError(f"Option '{option_id}' must be >= {self.minimum}")
		if self.maximum is not None and value > self.maximum:
			raise ValueError(f"Option '{option_id}' must be <= {self.maximum}")


class ToolDefinition(BaseModel):
	"""Definition of one configurable CLI tool."""

	name: str
	category: str
	options: dict[str, ToolOptionDefinition]


class ToolCatalog(BaseModel):
	"""Repository-controlled configurable tool catalog."""

	tools: dict[str, ToolDefinition]


class FastqMatchResult(BaseModel):
	"""Deterministic FASTQ matching output for one sample sheet."""

	uploaded_samples: list[UploadedSample]
	extra_fastqs: list[str] = Field(default_factory=list)


class RunPaths(BaseModel):
	"""Paths created for a staged run directory."""

	run_id: str
	root: Path
	input_dir: Path
	uploads_dir: Path
	results_dir: Path
	run_json: Path
	samples_csv: Path
	config_yaml: Path


class RunDefinition(BaseModel):
	"""Run-scoped data model for prepared analyses."""

	run_id: str | None = None
	uploaded_samples: list[UploadedSample]
	analyses: list[AnalysisDefinition]
	run_overrides: dict[str, dict[str, Any]] = Field(default_factory=dict)

	@model_validator(mode="after")
	def _validate_consistency(self) -> "RunDefinition":
		uploaded_samples = self.uploaded_samples or []
		analyses = self.analyses or []
		sample_ids = [sample.sample_id for sample in uploaded_samples]
		analysis_ids = [analysis.analysis_id for analysis in analyses]
		if len(sample_ids) != len(set(sample_ids)):
			raise ValueError("uploaded_samples must have unique sample_id values")
		if len(analysis_ids) != len(set(analysis_ids)):
			raise ValueError("analyses must have unique analysis_id values")
		unknown_samples = sorted(
			{analysis.uploaded_sample_id for analysis in analyses if analysis.uploaded_sample_id not in set(sample_ids)}
		)
		if unknown_samples:
			raise ValueError(f"analyses reference unknown uploaded_sample_id values: {unknown_samples}")
		return self


class PipelineConfigResolver:
	"""Resolve run-wide and per-analysis tool configuration for Snakemake."""

	def __init__(
		self,
		config: Mapping[str, Any],
		tool_catalog: ToolCatalog | None = None,
	):
		self.config = dict(config)
		self.tool_catalog = tool_catalog or load_tool_catalog()
		self.run_overrides = self._normalize_override_map(self.config.get("run_overrides") or {}, "run_overrides")
		analysis_overrides = self.config.get("analysis_overrides") or {}
		if not isinstance(analysis_overrides, Mapping):
			raise ValueError("analysis_overrides must be a mapping of analysis_id to overrides")
		self.analysis_overrides = {
			analysis_id: self._normalize_override_map(overrides, f"analysis_overrides.{analysis_id}")
			for analysis_id, overrides in analysis_overrides.items()
		}

	def get_analysis_parameter(self, analysis_id: str, tool_id: str, option_id: str) -> Any:
		"""Return the effective value for one analysis/tool option."""

		definition = self._get_option_definition(tool_id, option_id)
		analysis_values = self.analysis_overrides.get(analysis_id, {})
		if tool_id in analysis_values and option_id in analysis_values[tool_id]:
			return analysis_values[tool_id][option_id]
		if tool_id in self.run_overrides and option_id in self.run_overrides[tool_id]:
			return self.run_overrides[tool_id][option_id]
		for legacy_key in definition.legacy_keys:
			if legacy_key in self.config and not _is_missing(self.config[legacy_key]):
				return definition.validate_value(self.config[legacy_key], option_id)
		if definition.default is None:
			raise ValueError(f"No value found for option '{tool_id}.{option_id}'")
		return definition.validate_value(definition.default, option_id)

	def render_tool_args(
		self,
		analysis_id: str,
		tool_id: str,
		legacy_extra_key: str | None = None,
	) -> str:
		"""Render the effective CLI arguments for one tool."""

		if tool_id not in self.tool_catalog.tools:
			raise ValueError(f"Unknown tool_id '{tool_id}'")
		tokens: list[str] = []
		for option_id, definition in self.tool_catalog.tools[tool_id].options.items():
			value = self.get_analysis_parameter(analysis_id, tool_id, option_id)
			tokens.extend(definition.render_tokens(value))
		if legacy_extra_key:
			legacy_value = str(self.config.get(legacy_extra_key, "")).strip()
			if legacy_value:
				tokens.append(legacy_value)
		return " ".join(tokens)

	def _get_option_definition(self, tool_id: str, option_id: str) -> ToolOptionDefinition:
		if tool_id not in self.tool_catalog.tools:
			raise ValueError(f"Unknown tool_id '{tool_id}'")
		tool = self.tool_catalog.tools[tool_id]
		if option_id not in tool.options:
			raise ValueError(f"Unknown option '{tool_id}.{option_id}'")
		return tool.options[option_id]

	def _normalize_override_map(
		self,
		overrides: Mapping[str, Any],
		context: str,
	) -> dict[str, dict[str, Any]]:
		if not isinstance(overrides, Mapping):
			raise ValueError(f"{context} must be a mapping of tool IDs to option mappings")
		normalized: dict[str, dict[str, Any]] = {}
		for tool_id, tool_values in overrides.items():
			if tool_id not in self.tool_catalog.tools:
				raise ValueError(f"Unknown tool_id '{tool_id}' in {context}")
			if not isinstance(tool_values, Mapping):
				raise ValueError(f"{context}.{tool_id} must be a mapping of option IDs to values")
			normalized[tool_id] = {}
			for option_id, value in tool_values.items():
				definition = self._get_option_definition(tool_id, option_id)
				normalized[tool_id][option_id] = definition.validate_value(value, option_id)
		return normalized


def load_reference_catalog(path: Path | str = DEFAULT_REFERENCE_CATALOG_PATH) -> ReferenceCatalog:
	"""Load the repository reference catalog."""

	return _model_validate(ReferenceCatalog, yaml.safe_load(_normalize_path(path).read_text()) or {})


def load_tool_catalog(path: Path | str = DEFAULT_TOOL_CATALOG_PATH) -> ToolCatalog:
	"""Load the repository tool-option catalog."""

	return _model_validate(ToolCatalog, yaml.safe_load(_normalize_path(path).read_text()) or {})


def parse_illumina_samplesheet(path: Path | str) -> list[str]:
	"""Parse Illumina SampleSheet.csv sample IDs from the [Data] section."""

	text = _normalize_path(path).read_text()
	lines = text.splitlines()
	data_index = next((index for index, line in enumerate(lines) if line.strip().lower() == "[data]"), None)
	if data_index is not None:
		csv_text = "\n".join(line for line in lines[data_index + 1 :] if line.strip())
	else:
		csv_text = text
	try:
		table = pd.read_csv(StringIO(csv_text), dtype=str)
	except Exception as exc:
		raise ValueError(f"Unable to parse SampleSheet CSV: {exc}") from exc
	if table.empty:
		raise ValueError("SampleSheet contains no data rows")
	sample_column = next(
		(column for column in ("Sample_ID", "Sample_Name") if column in table.columns),
		None,
	)
	if sample_column is None:
		raise ValueError("SampleSheet must contain either Sample_ID or Sample_Name in the [Data] section")
	sample_ids = [str(value).strip() for value in table[sample_column].tolist() if not _is_missing(value)]
	if not sample_ids:
		raise ValueError("SampleSheet does not contain any usable sample IDs")
	duplicates = sorted({sample_id for sample_id in sample_ids if sample_ids.count(sample_id) > 1})
	if duplicates:
		raise ValueError(f"SampleSheet contains duplicate sample IDs: {duplicates}")
	return sample_ids


def match_fastqs_to_samples(
	sample_ids: Sequence[str],
	fastq_files: Sequence[Path | str],
) -> FastqMatchResult:
	"""Match FASTQ filenames to SampleSheet sample IDs using Illumina naming rules."""

	compiled_patterns = {
		sample_id: re.compile(
			rf"^{re.escape(sample_id)}_S\d+(?:_L\d{{3}})?_R(?P<read>[12])_001(?:\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)$"
		)
		for sample_id in sample_ids
	}
	per_sample: dict[str, dict[str, list[Path]]] = {
		sample_id: {"1": [], "2": []} for sample_id in sample_ids
	}
	seen_names: set[str] = set()
	extra_fastqs: list[str] = []
	for fastq_file in fastq_files:
		path = _normalize_path(fastq_file)
		if path.name in seen_names:
			raise ValueError(f"Duplicate FASTQ filename: '{path.name}'")
		seen_names.add(path.name)
		if FASTQ_NAME_PATTERN.fullmatch(path.name) is None:
			raise ValueError(f"Unsupported FASTQ filename: '{path.name}'")
		matches: list[tuple[str, str]] = []
		for sample_id, pattern in compiled_patterns.items():
			match = pattern.fullmatch(path.name)
			if match is not None:
				matches.append((sample_id, match.group("read")))
		if not matches:
			extra_fastqs.append(path.name)
			continue
		if len(matches) > 1:
			raise ValueError(f"Ambiguous FASTQ filename '{path.name}' matches multiple samples")
		sample_id, read = matches[0]
		per_sample[sample_id][read].append(path)
	uploaded_samples: list[UploadedSample] = []
	for sample_id, reads in per_sample.items():
		if not reads["1"]:
			raise ValueError(f"Missing R1 FASTQ for sample '{sample_id}'")
		if not reads["2"]:
			raise ValueError(f"Missing R2 FASTQ for sample '{sample_id}'")
		if len(reads["1"]) > 1:
			raise ValueError(f"Ambiguous R1 FASTQs for sample '{sample_id}'")
		if len(reads["2"]) > 1:
			raise ValueError(f"Ambiguous R2 FASTQs for sample '{sample_id}'")
		uploaded_samples.append(
			UploadedSample(sample_id=sample_id, r1=reads["1"][0], r2=reads["2"][0])
		)
	return FastqMatchResult(uploaded_samples=uploaded_samples, extra_fastqs=sorted(extra_fastqs))


def hash_files(file_names: Iterable[Path | str]) -> str:
	"""Hash one or more files using SHA256."""

	sha256 = hashlib.sha256()
	for file_name in file_names:
		with open(file_name, "rb") as handle:
			while True:
				data = handle.read(65536)
				if not data:
					break
				sha256.update(data)
	return sha256.hexdigest()


def build_pipeline_samples_map(sample_sheet: Path | str) -> dict[str, dict[str, Any]]:
	"""Load the current pipeline sample sheet into the Snakefile-compatible structure."""

	table = pd.read_csv(sample_sheet, sep=",", engine="python")
	required_columns = {"sample", "reference_fasta", "r1", "r2"}
	if not required_columns.issubset(set(table.columns)):
		raise ValueError("sample_sheet must contain columns: sample, reference_fasta, r1, r2")
	has_bed_column = "bed" in set(table.columns)
	samples: dict[str, dict[str, Any]] = {}
	for row in table.itertuples(index=False):
		sample = str(row.sample).strip()
		if _is_missing(sample):
			continue
		reference_fasta = str(row.reference_fasta).strip()
		r1 = str(row.r1).strip()
		r2 = str(row.r2).strip()
		bed = str(row.bed).strip() if has_bed_column else ""
		if _is_missing(reference_fasta) or _is_missing(r1) or _is_missing(r2):
			raise ValueError(f"Sample '{sample}' has missing reference_fasta, r1 or r2 in sample_sheet")
		if sample in samples:
			raise ValueError(f"Duplicate sample in sample_sheet: '{sample}'")
		reference_names = _read_fasta_names(reference_fasta)
		bed_names: list[str] = []
		if not _is_missing(bed):
			bed_names = _read_bed_names(bed)
			unknown_bed_names = sorted(set(bed_names) - set(reference_names))
			if unknown_bed_names:
				raise ValueError(
					f"BED file for sample '{sample}' contains references not present in FASTA: {unknown_bed_names}"
				)
		effective_references = bed_names if bed_names else reference_names
		if not effective_references:
			raise ValueError(f"Sample '{sample}' has no usable references after FASTA/BED parsing")
		files_to_hash = [reference_fasta] if _is_missing(bed) else [reference_fasta, bed]
		samples[sample] = {
			"reference_fasta": reference_fasta,
			"bed": "" if _is_missing(bed) else bed,
			"r1": r1,
			"r2": r2,
			"references": effective_references,
			"reference_hash": hash_files(files_to_hash),
		}
	if not samples:
		raise ValueError("sample_sheet was provided but no samples were parsed")
	return samples


def generate_samples_dataframe(
	analyses: Sequence[AnalysisDefinition],
	uploaded_samples: Sequence[UploadedSample],
	reference_catalog: ReferenceCatalog | None = None,
) -> pd.DataFrame:
	"""Create the legacy sample,reference_fasta,bed,r1,r2 CSV rows."""

	catalog = reference_catalog or load_reference_catalog()
	samples_by_id = {sample.sample_id: sample for sample in uploaded_samples}
	rows: list[dict[str, str]] = []
	for analysis in analyses:
		if analysis.uploaded_sample_id not in samples_by_id:
			raise ValueError(
				f"Analysis '{analysis.analysis_id}' references unknown uploaded sample '{analysis.uploaded_sample_id}'"
			)
		uploaded_sample = samples_by_id[analysis.uploaded_sample_id]
		reference_fasta, bed = catalog.resolve_paths(analysis.reference_id, analysis.primer)
		rows.append(
			{
				"sample": analysis.analysis_id,
				"reference_fasta": str(reference_fasta),
				"bed": bed,
				"r1": str(uploaded_sample.r1),
				"r2": str(uploaded_sample.r2),
			}
		)
	return pd.DataFrame(rows, columns=["sample", "reference_fasta", "bed", "r1", "r2"])


def write_samples_csv(
	output_path: Path | str,
	analyses: Sequence[AnalysisDefinition],
	uploaded_samples: Sequence[UploadedSample],
	reference_catalog: ReferenceCatalog | None = None,
) -> Path:
	"""Write the legacy pipeline sample CSV."""

	output = _normalize_path(output_path)
	dataframe = generate_samples_dataframe(analyses, uploaded_samples, reference_catalog)
	dataframe.to_csv(output, index=False)
	return output


def build_run_config(
	base_config: Mapping[str, Any],
	sample_sheet_path: Path | str,
	res_dir: Path | str,
	run_overrides: Mapping[str, Any] | None = None,
	analysis_overrides: Mapping[str, Any] | None = None,
	analyses: Sequence[AnalysisDefinition] | None = None,
	tool_catalog: ToolCatalog | None = None,
) -> dict[str, Any]:
	"""Build a run-specific config.yaml payload for Snakemake."""

	config = dict(base_config)
	config["sample_sheet"] = str(_normalize_path(sample_sheet_path))
	config["res_dir"] = str(_normalize_path(res_dir))
	config["run_overrides"] = dict(run_overrides or {})
	config["analysis_overrides"] = dict(analysis_overrides or _collect_analysis_overrides(analyses or []))
	PipelineConfigResolver(config, tool_catalog=tool_catalog)
	return config


def write_run_config(output_path: Path | str, config: Mapping[str, Any]) -> Path:
	"""Write a run-specific config.yaml file."""

	output = _normalize_path(output_path)
	output.write_text(yaml.safe_dump(dict(config), sort_keys=False))
	return output


def create_run_directory(
	storage_root: Path | str,
	run_id: str | None = None,
	created_at: datetime | None = None,
) -> RunPaths:
	"""Create the stable run directory used from staging through preparation."""

	root = _normalize_path(storage_root)
	effective_created_at = created_at or datetime.now(timezone.utc)
	timestamp = effective_created_at.strftime("%Y%m%d-%H%M%S")
	effective_run_id = run_id or f"{timestamp}-{token_hex(3)}"
	run_root = root / "runs" / effective_run_id
	input_dir = run_root / "input"
	uploads_dir = run_root / ".uploads"
	results_dir = run_root / "results"
	for directory in (input_dir, uploads_dir, results_dir):
		directory.mkdir(parents=True, exist_ok=False)
	run_json = run_root / "run.json"
	run_json.write_text(
		json.dumps(
			{
				"run_id": effective_run_id,
				"status": "staging",
				"created_at": effective_created_at.isoformat(),
			},
			indent=2,
		)
	)
	return RunPaths(
		run_id=effective_run_id,
		root=run_root,
		input_dir=input_dir,
		uploads_dir=uploads_dir,
		results_dir=results_dir,
		run_json=run_json,
		samples_csv=run_root / "samples.csv",
		config_yaml=run_root / "config.yaml",
	)


def generate_snakemake_command(
	repo_root: Path | str,
	config_path: Path | str,
	cores: int = 60,
	use_conda: bool = True,
	use_singularity: bool = True,
	singularity_bind_paths: Sequence[Path | str] | None = None,
	extra_args: Sequence[str] | None = None,
) -> str:
	"""Build the Snakemake command shown to the user for a prepared run."""

	command = [
		"snakemake",
		"--snakefile",
		str(_normalize_path(repo_root) / "Snakefile"),
		"--configfile",
		str(_normalize_path(config_path)),
		"--cores",
		str(cores),
		"--printshellcmds",
		"--keep-incomplete",
	]
	if use_conda:
		command.append("--use-conda")
	if use_singularity:
		command.append("--use-singularity")
		if singularity_bind_paths:
			joined_paths = ",".join(str(_normalize_path(path)) for path in singularity_bind_paths)
			command.extend(["--singularity-args", f"--bind {joined_paths}"])
	if extra_args:
		command.extend(extra_args)
	return " ".join(shlex.quote(part) for part in command)


def _read_fasta_names(fasta_path: Path | str) -> list[str]:
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


def _read_bed_names(bed_path: Path | str) -> list[str]:
	names: list[str] = []
	seen: set[str] = set()
	with open(bed_path) as handle:
		reader = csv.reader(handle, delimiter="\t")
		for row in reader:
			if not row:
				continue
			first = row[0].strip()
			if not first or first.startswith("#") or first in {"track", "browser"}:
				continue
			if first not in seen:
				seen.add(first)
				names.append(first)
	return names


def _collect_analysis_overrides(analyses: Sequence[AnalysisDefinition]) -> dict[str, dict[str, Any]]:
	collected: dict[str, dict[str, Any]] = {}
	for analysis in analyses:
		if analysis.tool_overrides:
			collected[analysis.analysis_id] = analysis.tool_overrides
	return collected