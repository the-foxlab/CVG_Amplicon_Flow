# GitHub Copilot Instructions

## Project goal

This repository contains an existing Snakemake bioinformatics pipeline.

We are adding a browser-based application that replaces the current workflow of manually uploading Illumina sequencing data to Galaxy.

The application runs on the same Linux workstation that runs Snakemake.

Users access it from a Windows computer on the institute's local network.

The typical workflow is:

```text
Illumina
   |
   v
USB drive
   |
   v
Windows computer
   |
   v
browser
   |
   v
Linux workstation
   |
   v
Snakemake
```

The MVP must allow the user to:

1. Upload an Illumina `SampleSheet.csv`.
2. Select FASTQ files.
3. Validate that the FASTQs correspond to the SampleSheet.
4. Upload FASTQs resumably to the Linux workstation.
5. Configure one or more analyses for each uploaded sample.
6. Select reference FASTA and primer scheme for each analysis.
7. Optionally configure arbitrary supported parameters of individual command-line tools on a per-analysis basis.
8. Optionally configure parameters for the whole run.
9. Generate the existing Snakemake `samples.csv`.
10. Generate a run-specific `config.yaml`.
11. Show the exact command required to run Snakemake.

The MVP does not execute Snakemake.

The architecture should make pipeline execution and monitoring straightforward to add later.

---

# Development priority

There is limited development time.

Prioritize a working end-to-end MVP over unnecessary abstraction.

Do not introduce features that are not required by this specification.

Do not perform unrelated refactoring.

Reuse the existing pipeline structure wherever practical.

Before changing the pipeline:

1. inspect the existing Snakefile and rule files
2. inspect the existing `config.yaml`
3. inspect existing sample sheets
4. inspect all existing command-line tool invocations
5. identify how configuration values are currently consumed

Treat the current pipeline behavior as authoritative unless this specification explicitly requires an extension.

Implement incrementally, but do not stop after producing only scaffolding or architectural proposals when implementation has been requested.

---

# General Python style

* Prefer simple and explicit Python.
* Keep functions focused and independently testable.
* Prefer ordinary functions and small data models over elaborate class hierarchies.
* Add abstractions only when they provide concrete value.
* Use type hints.
* Use `pathlib.Path` for filesystem paths.
* Use Python logging rather than `print()` for diagnostics.
* Raise meaningful exceptions.
* Avoid broad `except Exception` unless justified.
* Prefer standard-library functionality where practical.
* Avoid unnecessary third-party dependencies.
* Use Pydantic where typed validation materially simplifies the implementation.
* Use pytest for testing.

---

# Application technology

Use:

* Python
* FastAPI
* Pydantic
* Jinja2
* HTMX
* PyYAML
* pytest

For FASTQ upload use:

* Uppy
* `@uppy/core`
* `@uppy/dashboard`
* `@uppy/tus`
* tus resumable-upload protocol
* tusd v2
* nginx

Do not introduce:

* React
* Vue
* Angular
* a frontend build framework beyond what is necessary to package Uppy
* a relational database
* a job queue
* custom resumable-upload protocol
* custom FASTQ chunking implementation

Serve required frontend assets locally in production.

Do not depend on a public CDN at runtime.

Pin frontend dependency versions.

---

# Architectural boundaries

Keep these concerns separated:

```text
Web interface
    |
    +--> run management
    |
    +--> upload coordination
    |
    v
Domain logic
    |
    +--> SampleSheet parsing
    +--> FASTQ matching
    +--> sample models
    +--> analysis models
    +--> reference/primer catalog
    +--> tool option definitions
    +--> configuration resolution
    +--> validation
    +--> samples.csv generation
    +--> config.yaml generation
    +--> Snakemake command generation
```

Domain logic must not depend on FastAPI request objects, Jinja templates, Uppy, tus, or tusd.

Do not put substantial pipeline logic directly inside route handlers.

---

# Terminology

Distinguish carefully between an uploaded sample and an analysis.

## Uploaded sample

An uploaded sample originates from the Illumina SampleSheet and corresponds to an R1/R2 FASTQ pair.

Conceptually:

```python
UploadedSample(
    sample_id="1545554-HSV2_S82",
    r1=Path(...),
    r2=Path(...),
)
```

## Analysis

One uploaded sample can produce one or more analyses.

Each analysis corresponds to one row in the generated Snakemake sample sheet and therefore one value of the Snakemake `{sample}` wildcard.

For example:

```text
1545554-HSV2_S82
    |
    +--> 1545554-HSV2_S82_trimmed
    |
    +--> 1545554-HSV2_S82_not_trimmed
```

The FASTQ files must not be duplicated when multiple analyses use the same uploaded sample.

Use the term `analysis_id` internally for the value that will become the Snakemake `{sample}` wildcard.

---

# Existing Snakemake sample sheet

The current Snakemake sample-sheet schema is:

```csv
sample,reference_fasta,bed,r1,r2
```

Example:

```csv
1545554-HSV2_S82_trimmed,/path/refs/HSV_2_both.fasta,/path/beds/HSV_2_both.bed,/path/fastqs/1545554-HSV2_S82_L001_R1_001.fastq.gz,/path/fastqs/1545554-HSV2_S82_L001_R2_001.fastq.gz
1545554-HSV2_S82_not_trimmed,/path/refs/HSV_2_both.fasta,,/path/fastqs/1545554-HSV2_S82_L001_R1_001.fastq.gz,/path/fastqs/1545554-HSV2_S82_L001_R2_001.fastq.gz
```

Preserve compatibility with this schema for the MVP.

The GUI must not expose raw FASTA, BED, or FASTQ filesystem paths to normal users.

Internally, use stable identifiers and translate them to filesystem paths while generating `samples.csv`.

---

# Reference and primer catalog

Maintain a repository-controlled catalog describing references and compatible primer schemes.

A YAML representation is appropriate.

Conceptually:

```yaml
references:

  hsv2_both:
    name: HSV-2 both
    fasta: /path/to/refs/HSV_2_both.fasta

    primer_schemes:
      hsv2_both:
        name: HSV-2 both primers
        bed: /path/to/beds/HSV_2_both.bed
```

Each reference must have:

* stable identifier
* human-readable name
* FASTA path

Each primer scheme must have:

* stable identifier
* human-readable name
* BED path

Primer schemes belong to, or explicitly declare compatibility with, references.

The GUI must only allow valid reference and primer combinations.

Do not hard-code references or primer schemes into HTML or JavaScript.

---

# Analysis creation

For every uploaded Illumina sample, the user configures one or more analyses.

At minimum an analysis contains:

```text
analysis_id
uploaded sample
reference
primer processing
tool configuration overrides
```

Primer processing must explicitly represent either:

```text
trim with selected primer scheme
```

or:

```text
do not perform primer trimming
```

An empty BED path in the generated CSV represents no primer trimming.

If useful, the GUI may provide a convenience option to generate both:

```text
<sample>_trimmed
<sample>_not_trimmed
```

from the same uploaded sample.

Analysis identifiers should normally be generated by the application rather than manually entered.

They must be:

* unique within the run
* filesystem safe
* valid as Snakemake wildcard values

---

# Configuration model

There are three configuration levels:

```text
pipeline default
      |
      v
optional run-wide override
      |
      v
optional per-analysis override
```

The most specific explicitly configured value wins.

Conceptually:

```text
analysis override
    > run override
    > pipeline default
```

Per-analysis overrides must be keyed by `analysis_id`, because this is the Snakemake `{sample}` wildcard.

Example:

```yaml
analysis_overrides:

  1545554-HSV2_S82_trimmed:
    consensus:
      min_depth: 50

  1545554-HSV2_S82_not_trimmed:
    consensus:
      min_frequency: 0.7
```

Do not duplicate inherited defaults into every analysis.

---

# Tool configuration

The application must support configuration of command-line options for the tools used by the Snakemake pipeline.

This includes tool-specific options for mapping, trimming, consensus generation, variant calling, filtering, and other command-line tools currently used by the workflow.

The architecture must not be limited to a predefined small list such as depth or allele frequency.

Instead, configuration must be driven by structured definitions of individual tool options.

## Critical rule

Never allow the browser to supply an arbitrary command-line fragment such as:

```text
-k 3 -n 4 -c foo
```

Never expose a generic:

```text
extra arguments
```

free-text box.

Never directly interpolate arbitrary user-provided shell text into a Snakemake `shell` command.

Each configurable command-line option must be represented individually and validated according to its type.

Examples:

```text
minimum depth
    integer input

minimum allele frequency
    float input

mapping quality
    integer input

algorithm/mode
    select control

enable feature
    checkbox

reference preset
    select control
```

For a CLI option such as:

```text
-k INT
```

represent:

```text
option: -k
type: integer
value: 17
```

rather than storing:

```text
"-k 17"
```

For a boolean CLI flag such as:

```text
-Y
```

represent it as a boolean.

For a fixed-choice argument, use an enum/select.

For repeatable arguments, use a typed list where necessary.

Do not allow the user to configure:

* the executable path
* shell redirections
* pipes
* command separators
* environment-variable assignments
* arbitrary shell syntax

unless explicitly added in a future requirement.

---

# Tool option catalog

Create a repository-controlled structured catalog defining configurable CLI options.

Do not hard-code individual CLI options throughout route handlers and HTML templates.

The exact representation can be YAML or typed Python data, whichever results in the simpler implementation.

Conceptually:

```yaml
tools:

  bwa_mem:
    name: BWA MEM

    options:

      min_seed_length:
        flag: "-k"
        label: Minimum seed length
        type: int
        min: 1

      band_width:
        flag: "-w"
        label: Band width
        type: int
        min: 1

      mark_shorter_split_hits:
        flag: "-M"
        label: Mark shorter split hits
        type: bool

  ivar_consensus:
    name: iVar consensus

    options:

      min_depth:
        flag: "-m"
        label: Minimum depth
        type: int
        min: 0

      min_frequency:
        flag: "-t"
        label: Minimum consensus frequency
        type: float
        min: 0
        max: 1
```

This example is illustrative.

Inspect the actual pipeline commands and installed tool versions before defining real options.

Tool definitions should support at least:

```text
stable option identifier
human-readable label
CLI flag
type
default/inherited value
optional minimum
optional maximum
optional enum choices
optional help text
optional repeatability
```

Supported types should include where useful:

```text
bool
int
float
string
enum
list
```

String-valued options are values for one known option.

They are not arbitrary command fragments.

---

# Scope of tool options

The GUI architecture must be capable of exposing all meaningful configurable options of each supported command-line tool individually.

Do not build a hard-coded special case only for:

```text
AF
depth
MAPQ
```

These are simply examples.

For every command-line tool currently invoked by the Snakemake workflow:

1. identify the command/version
2. identify its configurable CLI options
3. represent supported options in the tool option catalog
4. allow those options to be set at run or analysis level where meaningful

It is acceptable for the initial implementation to define these catalogs manually.

Do not attempt to dynamically parse command-line `--help` output at runtime.

The catalog is repository-controlled and versioned together with the pipeline.

---

# Tool option GUI

For every tool, show an advanced configuration section.

Conceptually:

```text
Advanced settings

Mapping
  BWA MEM

    Minimum seed length
    inherited: 19
    [ ] Override
        [ 17 ]

    Band width
    inherited: 100
    [ ] Override
        [ 80 ]

    Mark shorter split hits
    inherited: off
    [ ] Override
        [x]

Consensus
  iVar consensus

    Minimum depth
    inherited: 20
    [ ] Override
        [ 50 ]

    Minimum frequency
    inherited: 0.8
    [ ] Override
        [ 0.9 ]
```

The GUI must clearly distinguish:

```text
inherited value
```

from:

```text
explicit override
```

If an override is disabled, do not write a redundant override to the generated configuration.

Tool settings should normally be collapsed under an "Advanced settings" section.

The normal workflow should remain simple.

---

# Run-level versus analysis-level tool configuration

Allow applicable tool settings at both:

```text
run level
analysis level
```

Run-level changes apply to every analysis unless overridden.

Analysis-level settings apply only to one `analysis_id`.

Example:

```text
pipeline default depth: 20

run override: 30

analysis A: inherited -> 30
analysis B: override 50 -> 50
```

The effective configuration resolver must be centralized and tested.

Do not scatter nested `config.get()` logic throughout rules.

---

# Snakemake configuration integration

Extend the existing Snakemake pipeline so rules can resolve configuration by wildcard.

Use wildcard-dependent `params` functions where required.

Conceptually:

```python
params:
    min_depth=lambda wildcards: get_analysis_parameter(
        wildcards.sample,
        "ivar_consensus",
        "min_depth",
    )
```

Create one centralized helper that resolves:

```text
analysis override
    |
    v
run override
    |
    v
pipeline default
```

Do not duplicate resolution logic in each rule.

Keep current behavior unchanged when no overrides are specified.

Existing runs without `analysis_overrides` must continue to work if practical.

---

# Safe command construction

All tool command arguments controlled by the GUI must originate from validated structured configuration.

Do not build a shell argument by concatenating arbitrary browser-provided strings.

When possible:

* render known flags from trusted tool definitions
* validate typed values
* quote values safely when shell invocation requires it

The browser may provide:

```json
{
  "min_seed_length": 17
}
```

The trusted server-side catalog determines that this corresponds to:

```text
-k 17
```

The browser must not determine that `-k` is the associated CLI flag.

---

# File upload architecture

FASTQ upload is part of the MVP.

FASTQ files may be many gigabytes.

Do not upload them using an ordinary custom FastAPI multipart implementation.

Use:

```text
Browser:
    Uppy
        @uppy/core
        @uppy/dashboard
        @uppy/tus

Protocol:
    tus

Upload server:
    tusd v2

Application:
    FastAPI

Reverse proxy:
    nginx

Storage:
    Linux filesystem
```

FastAPI is the control plane.

tusd is the large-file data plane.

FASTQ bytes must not be proxied through Python.

---

# Reverse proxy

Conceptually:

```text
Browser
   |
   v
nginx
   |
   +--> /files/ --> tusd
   |
   +--> /       --> FastAPI
```

FastAPI and tusd should normally listen only on localhost.

nginx is the service exposed to the institute network.

For `/files/`:

* disable request buffering
* permit large request bodies
* proxy directly to tusd

The internal tusd hook endpoint must not be exposed publicly.

---

# Run lifecycle

Create a run before FASTQ upload.

At minimum:

```text
staging
prepared
```

A run becomes `prepared` only after:

* SampleSheet parsing succeeds
* upload plan validation succeeds
* required FASTQs are fully uploaded
* FASTQ matching succeeds
* analyses are valid
* reference/primer selections are valid
* tool configurations are valid
* `samples.csv` is generated
* `config.yaml` is generated

The MVP does not execute Snakemake.

---

# Run directory

Use one stable run directory from upload through finalization.

Conceptually:

```text
<run-storage>/
    runs/
        <run-id>/

            input/
                SampleSheet.csv
                SAMPLE_R1.fastq.gz
                SAMPLE_R2.fastq.gz

            .uploads/
                upload-plan.json
                <upload-id>.json

            run.json
            samples.csv
            config.yaml

            results/
```

Avoid moving or copying multi-gigabyte FASTQs during finalization.

The storage root must be configurable.

Do not store sequencing data inside the Git repository.

---

# SampleSheet upload

The SampleSheet is small.

Upload it through a normal FastAPI `UploadFile`.

Parse it before starting FASTQ transfer.

Reject malformed SampleSheets before FASTQ upload.

---

# FASTQ pre-flight validation

Uppy must not automatically start uploading immediately after file selection.

Before binary transfer, send only:

```text
filename
size
```

to FastAPI.

FastAPI must validate the planned upload.

At minimum validate:

* accepted FASTQ suffix
* no duplicate filenames
* no zero-byte files
* expected R1 files
* expected R2 files
* missing pairs
* ambiguous matching
* extra files
* total expected size
* available storage space where practical

Do not transfer gigabytes of data when the selected set is already known to be invalid.

---

# Uppy configuration

Use Uppy Dashboard for:

* multi-file selection
* progress
* pause
* resume
* retry
* failure display

Use `@uppy/tus` for transfer.

Initial concurrent upload limit:

```text
3
```

Do not define a custom tus chunk size unless testing shows a concrete requirement.

Automatic retry must remain enabled.

Associate uploads with trusted metadata including:

```text
run_id
filename
```

---

# tusd hooks

Integrate tusd with FastAPI via HTTP hooks.

At minimum use:

```text
pre-create
post-finish
```

## pre-create

Validate:

* run exists
* state is `staging`
* filename is present
* filename is safe
* filename is present in the approved upload plan
* declared size matches the upload plan
* destination is unique
* the file has not already completed

Reject:

* absolute paths
* `..`
* `/`
* `\`
* empty names
* duplicate destination names

Never use an unvalidated browser filename directly as a path.

Use generated tus upload IDs.

Where appropriate, use the tusd hook mechanism to place the completed FASTQ directly in:

```text
<run>/input/<validated-original-filename>
```

## post-finish

Treat this as the authoritative signal that a file completed.

Record:

* upload ID
* run ID
* filename
* final path
* expected size
* completed size
* completion time

The operation must be idempotent.

---

# Upload completion

Do not consider an upload complete merely because:

* Uppy displays 100 percent
* a file exists
* its size appears correct

Use server-confirmed tus completion recorded by `post-finish`.

Finalization must remain blocked until all required files are complete.

---

# FASTQ matching

FASTQ matching must remain independent from Uppy and tus.

Match original validated filenames against SampleSheet entries using deterministic Illumina naming logic.

Do not use arbitrary substring matching if a deterministic rule is available.

Identify:

* matched R1
* matched R2
* missing files
* ambiguous files
* duplicate matches
* unmatched extra uploads

---

# Testing

Use small fixtures only.

Never add real large sequencing files to the test suite.

At minimum test:

## SampleSheet

* valid SampleSheet
* multiple samples
* malformed SampleSheet

## FASTQ matching

* valid R1/R2 pair
* missing R1
* missing R2
* ambiguous match
* extra FASTQ
* duplicate filename

## Upload control

* valid upload plan
* invalid upload plan
* unsafe filename
* wrong size
* unknown run
* `pre-create`
* `post-finish`
* duplicate `post-finish`
* incomplete upload prevents finalization

## References and primers

* valid reference
* valid primer scheme
* invalid combination
* no-trimming analysis

## Analyses

* one uploaded sample produces one analysis
* one uploaded sample produces multiple analyses
* FASTQ paths are reused rather than duplicated
* analysis IDs are unique

## Tool options

* typed integer option
* typed float option
* boolean flag
* enum option
* invalid range
* unsupported option
* arbitrary option key rejection
* arbitrary raw CLI text rejection
* run-level inheritance
* analysis-level override
* analysis-level override wins over run setting
* run setting wins over pipeline default

## Output

* correct `samples.csv`
* correct `config.yaml`
* correct analysis override structure
* correct Snakemake command

---

# Working procedure for Copilot

For each substantial implementation task:

1. inspect relevant existing files
2. identify current behavior
3. implement the smallest coherent end-to-end change
4. add or update tests
5. run tests
6. report changed files and assumptions

Do not spend substantial time proposing alternative architectures after this specification has established the architecture.

When an implementation detail is unspecified, prefer the simplest solution consistent with this document and the existing codebase.

Do not create custom agent infrastructure, databases, frontend frameworks, or other architectural components unless required to make the specified MVP work.
