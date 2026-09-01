const galaxyPortRun = {
	uppy: null,
	planOk: false,
	analysisCounter: 0,

	init() {
		if (window.galaxyPortRunId) {
			this.initUppy();
			this.ensureInitialAnalysis();
		}
	},

	initUppy() {
		if (!window.Uppy || !document.getElementById("uppy-dashboard")) {
			return;
		}
		this.uppy = new Uppy.Core({
			autoProceed: false,
			restrictions: {
				minNumberOfFiles: 1,
			},
		});
		this.uppy.use(Uppy.Dashboard, {
			target: "#uppy-dashboard",
			inline: true,
			height: 350,
			hideUploadButton: true,
			note: "Select FASTQs, validate plan, then start upload.",
		});
		this.uppy.use(Uppy.Tus, {
			endpoint: "/files/",
			limit: 3,
			retryDelays: [0, 1000, 3000, 5000],
		});
		this.uppy.on("file-added", (file) => {
			this.uppy.setFileMeta(file.id, {
				run_id: window.galaxyPortRunId,
				filename: file.name,
			});
		});
	},

	collectSelectedFiles() {
		if (!this.uppy) {
			return [];
		}
		return this.uppy.getFiles().map((file) => ({
			name: file.name,
			size: file.size,
		}));
	},

	async validatePlan(runId) {
		const selectedFiles = this.collectSelectedFiles();
		const response = await fetch(`/runs/${runId}/upload-plan`, {
			method: "POST",
			headers: {"Content-Type": "application/json"},
			body: JSON.stringify({files: selectedFiles}),
		});
		const output = document.getElementById("plan-output");
		const startButton = document.getElementById("start-upload-button");
		if (!response.ok) {
			const error = await response.json();
			output.textContent = `Validation failed: ${error.detail}`;
			this.planOk = false;
			if (startButton) {
				startButton.disabled = true;
			}
			return;
		}
		const plan = await response.json();
		window.galaxyPortPlan = plan;
		output.textContent = JSON.stringify(plan, null, 2);
		this.planOk = true;
		if (startButton) {
			startButton.disabled = false;
		}
	},

	startUpload() {
		if (!this.uppy || !this.planOk) {
			return;
		}
		this.uppy.upload();
	},

	ensureInitialAnalysis() {
		if (!document.getElementById("analysis-rows")) {
			return;
		}
		if ((window.galaxyPortSampleIds || []).length === 0) {
			return;
		}
		if (document.querySelectorAll(".analysis-row").length === 0) {
			this.addAnalysisRowDom(this.analysisCounter);
			this.analysisCounter += 1;
		}
	},

	addAnalysisRowDom(index) {
		const container = document.getElementById("analysis-rows");
		if (!container) {
			return;
		}
		const references = window.galaxyPortReferenceOptions || [];
		if (references.length === 0) {
			return;
		}
		const sampleId = (window.galaxyPortSampleIds || ["sample"])[0];
		const ref = references[0];
		const primerOptions = (ref.primer_schemes || []).map((primer) => `<option value="${primer.primer_scheme_id}">${primer.name}</option>`).join("");
		const html = `
		<div class="analysis-row" data-analysis-index="${index}">
			<h4>Analysis ${index + 1}</h4>
			<label>Analysis ID
				<input type="text" class="analysis-id" value="${sampleId}_${index + 1}" required>
			</label>
			<label>Uploaded Sample
				<select class="uploaded-sample-id">
					${(window.galaxyPortSampleIds || []).map((id) => `<option value="${id}">${id}</option>`).join("")}
				</select>
			</label>
			<label>Reference
				<select class="reference-id" onchange="galaxyPortRun.updatePrimerOptions(this)">
					${references.map((entry) => `<option value="${entry.reference_id}">${entry.name}</option>`).join("")}
				</select>
			</label>
			<label>Primer Mode
				<select class="primer-mode" onchange="galaxyPortRun.togglePrimerMode(this)">
					<option value="trim">Trim primers</option>
					<option value="none">Do not trim primers</option>
				</select>
			</label>
			<label class="primer-scheme-wrap">Primer Scheme
				<select class="primer-scheme-id">${primerOptions}</select>
			</label>
			<details>
				<summary>Per-analysis advanced settings</summary>
				${this.renderToolControls()}
			</details>
		</div>`;
		container.insertAdjacentHTML("beforeend", html);
	},

	renderToolControls() {
		return (window.galaxyPortToolCatalog || []).map((tool) => {
			const options = (tool.options || []).map((option) => {
				let valueControl = "";
				if (option.type === "bool") {
					valueControl = `<label>Value<input type="checkbox" class="option-value-bool" ${option.default ? "checked" : ""}></label>`;
				} else if (option.type === "enum") {
					const enumChoices = (option.choices || []).map((choice) => `<option value="${choice}" ${choice === option.default ? "selected" : ""}>${choice}</option>`).join("");
					valueControl = `<select class="option-value-enum">${enumChoices}</select>`;
				} else {
					valueControl = `<input type="text" class="option-value-text" value="${option.default}">`;
				}
				return `<div class="option-row" data-option-id="${option.option_id}" data-option-type="${option.type}">
					<label><input type="checkbox" class="option-enabled">Override ${option.label} (default: ${option.default})</label>
					${valueControl}
				</div>`;
			}).join("");
			return `<div class="tool-block" data-tool-id="${tool.tool_id}"><h5>${tool.category}: ${tool.name}</h5>${options}</div>`;
		}).join("");
	},

	updatePrimerOptions(referenceSelect) {
		const row = referenceSelect.closest(".analysis-row");
		if (!row) {
			return;
		}
		const primerSelect = row.querySelector(".primer-scheme-id");
		if (!primerSelect) {
			return;
		}
		const referenceId = referenceSelect.value;
		const reference = (window.galaxyPortReferenceOptions || []).find((entry) => entry.reference_id === referenceId);
		const primerOptions = (reference?.primer_schemes || []).map((primer) => `<option value="${primer.primer_scheme_id}">${primer.name}</option>`).join("");
		primerSelect.innerHTML = primerOptions;
	},

	togglePrimerMode(primerModeSelect) {
		const row = primerModeSelect.closest(".analysis-row");
		if (!row) {
			return;
		}
		const wrap = row.querySelector(".primer-scheme-wrap");
		if (!wrap) {
			return;
		}
		wrap.style.display = primerModeSelect.value === "trim" ? "grid" : "none";
	},

	parseTypedOption(optionRow) {
		const type = optionRow.dataset.optionType;
		if (type === "bool") {
			return !!optionRow.querySelector(".option-value-bool")?.checked;
		}
		if (type === "enum") {
			return optionRow.querySelector(".option-value-enum")?.value;
		}
		const textValue = optionRow.querySelector(".option-value-text")?.value ?? "";
		if (type === "int") {
			return Number.parseInt(textValue, 10);
		}
		if (type === "float") {
			return Number.parseFloat(textValue);
		}
		if (type === "list") {
			return textValue.split(",").map((item) => item.trim()).filter(Boolean);
		}
		return textValue;
	},

	collectOverrides(container) {
		const overrides = {};
		container.querySelectorAll(".tool-block").forEach((toolBlock) => {
			const toolId = toolBlock.dataset.toolId;
			toolBlock.querySelectorAll(".option-row").forEach((optionRow) => {
				const enabled = optionRow.querySelector(".option-enabled")?.checked;
				if (!enabled) {
					return;
				}
				const optionId = optionRow.dataset.optionId;
				if (!overrides[toolId]) {
					overrides[toolId] = {};
				}
				overrides[toolId][optionId] = this.parseTypedOption(optionRow);
			});
		});
		return overrides;
	},

	collectAnalyses() {
		const analyses = [];
		document.querySelectorAll(".analysis-row").forEach((row) => {
			const primerMode = row.querySelector(".primer-mode")?.value || "none";
			analyses.push({
				analysis_id: row.querySelector(".analysis-id")?.value,
				uploaded_sample_id: row.querySelector(".uploaded-sample-id")?.value,
				reference_id: row.querySelector(".reference-id")?.value,
				primer_mode: primerMode,
				primer_scheme_id: primerMode === "trim" ? row.querySelector(".primer-scheme-id")?.value : null,
				tool_overrides: this.collectOverrides(row),
			});
		});
		return analyses;
	},

	collectRunOverrides() {
		const wrap = document.getElementById("run-overrides");
		if (!wrap) {
			return {};
		}
		return this.collectOverrides(wrap);
	},

	async prepareRun(runId) {
		const payload = {
			analyses: this.collectAnalyses(),
			run_overrides: this.collectRunOverrides(),
		};
		const response = await fetch(`/runs/${runId}/prepare`, {
			method: "POST",
			headers: {"Content-Type": "application/json"},
			body: JSON.stringify(payload),
		});
		const output = document.getElementById("prepare-output");
		if (!response.ok) {
			const error = await response.json();
			output.textContent = `Preparation failed: ${error.detail}`;
			return;
		}
		const result = await response.json();
		output.textContent = JSON.stringify(result, null, 2);
		window.location.href = `/runs/${runId}`;
	},
};

window.galaxyPortRun = galaxyPortRun;

document.addEventListener("DOMContentLoaded", () => {
	galaxyPortRun.init();
	document.body.addEventListener("htmx:afterSwap", () => {
		if (document.querySelectorAll(".analysis-row").length > galaxyPortRun.analysisCounter) {
			galaxyPortRun.analysisCounter = document.querySelectorAll(".analysis-row").length;
		}
	});
	document.body.addEventListener("click", (event) => {
		const target = event.target;
		if (!(target instanceof HTMLElement)) {
			return;
		}
		if (target.getAttribute("hx-get")?.includes("analysis-row")) {
			const nextIndex = document.querySelectorAll(".analysis-row").length;
			target.setAttribute("hx-get", `/runs/${window.galaxyPortRunId}/analysis-row?index=${nextIndex}`);
		}
	});
});
