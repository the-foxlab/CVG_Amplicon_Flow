# Snakemake Run Preparation GUI

## MVP specification

## 1. Purpose

Build a browser-based interface for preparing sequencing analyses for the existing Snakemake pipeline.

The current workflow is:

```text
Illumina
   |
   v
USB drive
   |
   v
Windows workstation
   |
   v
upload to Galaxy
   |
   v
run analysis
   |
   v
download results
```

The new MVP should provide:

```text
Illumina
   |
   v
USB drive
   |
   v
Windows workstation
   |
   v
browser
   |
   v
Linux analysis workstation
   |
   +--> FASTQ storage
   +--> Snakemake samples.csv
   +--> Snakemake config.yaml
   |
   v
show Snakemake command
```

The MVP prepares the run.

It does not execute Snakemake yet.

---

# 2. Technology

Use:

```text
Backend:
    Python
    FastAPI
    Pydantic
    PyYAML

Frontend:
    Jinja2
    HTMX

Large-file upload:
    Uppy
        @uppy/core
        @uppy/dashboard
        @uppy/tus

Upload protocol:
    tus

Upload server:
    tusd v2

Reverse proxy:
    nginx

Testing:
    pytest
```

Do not build a separate SPA.

---

# 3. User workflow

The intended user workflow is:

```text
Create analysis run
        |
        v
Upload SampleSheet.csv
        |
        v
Parse Illumina samples
        |
        v
Select FASTQ files
        |
        v
Pre-flight filename validation
        |
        v
Upload FASTQs resumably
        |
        v
Verify server-side completion
        |
        v
Configure analyses
        |
        +--> reference
        +--> primer scheme
        +--> trimming / no trimming
        +--> per-analysis tool settings
        |
        v
Optional run-wide advanced settings
        |
        v
Validate
        |
        v
Generate samples.csv
        |
        v
Generate config.yaml
        |
        v
Display Snakemake command
```

The normal workflow should be simple.

Advanced configuration should be available but not required.

---

# 4. Uploaded sample versus analysis

These are distinct concepts.

## Uploaded sample

An Illumina sample corresponds to uploaded sequencing data:

```text
1545554-HSV2_S82

R1
1545554-HSV2_S82_L001_R1_001.fastq.gz

R2
1545554-HSV2_S82_L001_R2_001.fastq.gz
```

## Analysis

One uploaded sample may produce multiple Snakemake analyses.

Example:

```text
1545554-HSV2_S82
    |
    +--> 1545554-HSV2_S82_trimmed
    |
    +--> 1545554-HSV2_S82_not_trimmed
```

Both analyses use the same FASTQ files.

FASTQs must not be copied or uploaded twice.

Each analysis corresponds to one row in `samples.csv` and one `{sample}` wildcard value in Snakemake.

---

# 5. Existing sample-sheet format

The current pipeline expects:

```csv
sample,reference_fasta,bed,r1,r2
```

Example:

```csv
1545554-HSV2_S82_trimmed,/home/udo/shared/researchers/udo_gieraths/galaxy_port_data/refs/HSV_2_both.fasta,/home/udo/shared/researchers/udo_gieraths/galaxy_port_data/beds/HSV_2_both.bed,/home/udo/shared/researchers/udo_gieraths/galaxy_port_data/fastqs/1545554-HSV2_S82_L001_R1_001.fastq.gz,/home/udo/shared/researchers/udo_gieraths/galaxy_port_data/fastqs/1545554-HSV2_S82_L001_R2_001.fastq.gz

1545554-HSV2_S82_not_trimmed,/home/udo/shared/researchers/udo_gieraths/galaxy_port_data/refs/HSV_2_both.fasta,,/home/udo/shared/researchers/udo_gieraths/galaxy_port_data/fastqs/1545554-HSV2_S82_L001_R1_001.fastq.gz,/home/udo/shared/researchers/udo_gieraths/galaxy_port_data/fastqs/1545554-HSV2_S82_L001_R2_001.fastq.gz
```

Preserve this format.

Do not expose these filesystem paths directly in the GUI.

---

# 6. Reference catalog

Create a repository-controlled reference catalog.

For example:

```yaml
references:

  hsv2_both:
    name: HSV-2 both
    fasta: /home/udo/shared/researchers/udo_gieraths/galaxy_port_data/refs/HSV_2_both.fasta

    primer_schemes:

      hsv2_both:
        name: HSV-2 both primers
        bed: /home/udo/shared/researchers/udo_gieraths/galaxy_port_data/beds/HSV_2_both.bed
```

The user sees:

```text
HSV-2 both
```

rather than:

```text
/home/udo/.../HSV_2_both.fasta
```

References and primers must use stable IDs internally.

---

# 7. Analysis setup UI

After uploads complete, show the uploaded samples.

Conceptually:

```text
Sample                 R1              R2

1545554-HSV2_S82       uploaded ✓      uploaded ✓
```

Allow the user to create analysis entries.

For each analysis:

```text
Analysis ID
1545554-HSV2_S82_trimmed

Reference
[ HSV-2 both ▼ ]

Primer processing

(o) Trim primers
    Primer scheme:
    [ HSV-2 both primers ▼ ]

( ) Do not trim primers

Advanced settings
[ expand ]
```

Provide a convenience mechanism for generating both trimmed and untrimmed analyses from the same uploaded sample if practical.

---

# 8. Run-wide and analysis-specific configuration

Pipeline configuration has three levels:

```text
pipeline defaults
        |
        v
run-wide overrides
        |
        v
analysis-specific overrides
```

The most specific configured value wins.

Example:

```text
default minimum depth: 20

run override: 30

analysis A:
    no override
    effective value = 30

analysis B:
    override = 50
    effective value = 50
```

---

# 9. Tool configuration model

Users must be able to configure command-line tools used by the pipeline.

This should not be limited to a few parameters such as:

```text
allele frequency
depth
MAPQ
```

The system should support the different configurable options of each pipeline tool.

However, every CLI option must be represented individually.

Never allow:

```text
Additional arguments:
[ -k 3 -n 4 -c foo ]
```

Instead expose:

```text
Minimum seed length
[ 3 ]

Number of ...
[ 4 ]

Some mode
[ choice ▼ ]
```

The backend maps these typed values to trusted command-line flags.

---

# 10. Tool option definitions

Maintain a repository-controlled catalog of configurable options.

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

These examples are illustrative.

The actual catalog must be based on the tools and versions currently used by the Snakemake pipeline.

---

# 11. Supported option types

The tool-option system should support at least:

```text
bool
int
float
string
enum
list
```

Examples:

## Boolean flag

CLI:

```text
-M
```

GUI:

```text
[x] Mark shorter split hits
```

## Integer

CLI:

```text
-k 17
```

GUI:

```text
Minimum seed length
[ 17 ]
```

## Float

CLI:

```text
-t 0.8
```

GUI:

```text
Consensus frequency
[ 0.8 ]
```

## Enum

CLI accepts one of a predefined set.

GUI:

```text
Mode
[ mode A ▼ ]
```

## String

A known option may accept a string value.

This is allowed.

A string is the value of one known option.

It is not an arbitrary shell fragment.

---

# 12. Tool configuration safety

The browser must never determine the command syntax.

For example, the browser may send:

```json
{
    "min_seed_length": 17
}
```

The server-side trusted catalog determines that:

```text
min_seed_length
```

corresponds to:

```text
-k
```

and produces:

```text
-k 17
```

The user must not be able to provide:

```text
-k 17
```

directly.

The user must not be able to provide:

```text
-k 17 ; rm ...
```

or any arbitrary shell expression.

Do not expose:

```text
extra_args
raw_args
command
shell_options
```

free-text fields.

---

# 13. Tool coverage

Inspect every command-line tool invoked by the current workflow.

For every relevant tool:

1. identify the tool and version
2. identify command-line options that may reasonably be configured
3. create structured option definitions
4. allow these options to be configured in the advanced GUI
5. validate values according to their declared type and constraints

The architecture should permit essentially any CLI option of a supported tool to be represented without changing the application architecture.

Adding another option should normally require only adding another option definition to the catalog.

Do not dynamically parse CLI `--help` during production runtime.

The catalog is version-controlled with the pipeline.

---

# 14. Advanced configuration UI

The normal UI should not be overwhelmed with tool settings.

Use collapsible sections.

Conceptually:

```text
Advanced analysis settings

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


Consensus
  iVar

    Minimum depth
    inherited: 20
    [ ] Override
        [ 50 ]

    Consensus frequency
    inherited: 0.8
    [ ] Override
        [ 0.9 ]
```

Every field should show the inherited/effective value.

The user explicitly enables an override.

If an override is not enabled, do not store the value redundantly.

---

# 15. Run-wide advanced settings

Provide the same structured option system at run level.

Conceptually:

```text
Run-wide advanced settings

Mapping
    ...

Consensus
    ...

Variant calling
    ...
```

These become defaults for all analyses in the run.

Analysis-specific overrides take precedence.

---

# 16. Generated config.yaml

Represent overrides structurally.

Conceptually:

```yaml
run_overrides:

  ivar_consensus:
    min_depth: 30

analysis_overrides:

  1545554-HSV2_S82_trimmed:

    ivar_consensus:
      min_depth: 50
      min_frequency: 0.9

    bwa_mem:
      min_seed_length: 17

  1545554-HSV2_S82_not_trimmed:

    ivar_consensus:
      min_frequency: 0.7
```

The exact naming should be adapted to the existing pipeline where doing so avoids unnecessary changes.

Do not copy every pipeline default into this section.

Store explicit overrides only.

---

# 17. Snakemake integration

Extend the existing workflow with one central configuration resolver.

Conceptually:

```python
get_analysis_parameter(
    analysis_id,
    tool,
    option,
)
```

Resolution order:

```text
analysis override
        |
        v
run override
        |
        v
existing pipeline default
```

Rules should use Snakemake wildcard-aware `params` functions.

Conceptually:

```python
params:
    min_depth=lambda wildcards: get_analysis_parameter(
        wildcards.sample,
        "ivar_consensus",
        "min_depth",
    )
```

Do not duplicate this lookup logic across rules.

The existing pipeline should behave unchanged when no overrides are provided.

---

# 18. FASTQ upload

FASTQ files are uploaded directly from the Windows browser to the Linux workstation.

FASTQ files may be large.

Uploads must be resumable.

Use:

```text
Uppy
   |
   v
tus
   |
   v
nginx
   |
   v
tusd
   |
   v
Linux filesystem
```

FastAPI coordinates the run but does not receive FASTQ bytes.

---

# 19. Run creation

Create a run before upload.

States:

```text
staging
prepared
```

Example:

```text
/mnt/analysis/pipeline-runs/
    runs/
        20260831-122541-a1b2c3/
```

Store:

```text
input/
.uploads/
run.json
samples.csv
config.yaml
results/
```

FASTQs belong directly to this run.

Do not create another run directory and copy them during finalization.

---

# 20. SampleSheet upload

Upload `SampleSheet.csv` normally through FastAPI.

Store:

```text
<run>/input/SampleSheet.csv
```

Parse immediately.

If parsing fails, do not start FASTQ transfer.

---

# 21. FASTQ selection pre-flight

Use Uppy to select multiple `.fastq.gz` files.

Do not immediately start transfer.

First send only:

```text
filename
size
```

to FastAPI.

Perform matching before upload.

Show something like:

```text
Sample             R1        R2

ABC001             found     found
ABC002             found     missing
```

Show unmatched extra files separately.

Blocking errors must be fixed before binary upload starts.

---

# 22. Pre-flight validation

At minimum check:

```text
valid filename suffix
non-zero file size
no duplicate filenames
R1 present
R2 present where expected
no ambiguous matches
SampleSheet sample exists
available storage space where practical
```

Extra unrelated FASTQs may be warnings rather than errors.

---

# 23. Uppy behavior

Use Uppy Dashboard.

Provide:

```text
multi-file selection
per-file progress
overall progress
pause
resume
retry
failure information
```

Use `@uppy/tus`.

Initial concurrency:

```text
3 uploads
```

Keep automatic retries enabled.

Do not configure a custom chunk size unless testing provides a reason.

---

# 24. nginx and tusd

Browser-facing routing:

```text
/
    -> FastAPI

/files/
    -> tusd
```

For the tus upload route:

```text
disable request buffering
allow large request bodies
```

FastAPI and tusd should normally only listen on localhost.

The tusd hook endpoint must not be externally exposed.

---

# 25. tusd pre-create hook

Before accepting an upload validate:

```text
run exists
run is staging
filename is safe
filename is planned
size matches upload plan
destination is unique
upload not already completed
```

Reject filenames containing:

```text
..
/
\
absolute paths
```

Do not trust browser paths.

Preserve valid original FASTQ filenames.

---

# 26. tusd post-finish hook

`tusd` notifies FastAPI when a transfer completes.

This is the authoritative completion signal.

Persist:

```text
upload ID
run ID
filename
path
size
completion timestamp
```

The hook must be idempotent.

The run cannot be finalized until every required FASTQ is server-confirmed complete.

---

# 27. SampleSheet to FASTQ matching

Match FASTQs to Illumina samples using deterministic Illumina naming conventions.

Keep this code independent from tus.

The matching logic takes sample information plus filenames and returns:

```text
matched R1
matched R2
missing
ambiguous
unmatched extras
```

---

# 28. Final validation

A run may only become `prepared` when:

```text
SampleSheet valid
all required FASTQs uploaded
FASTQ matching valid
at least one valid analysis exists
references valid
primer combinations valid
analysis IDs unique
tool configuration valid
run-wide configuration valid
samples.csv successfully generated
config.yaml successfully generated
```

---

# 29. Generated samples.csv

For each analysis write:

```text
sample
reference_fasta
bed
r1
r2
```

For primer-trimmed analysis:

```text
bed = selected BED path
```

For untrimmed analysis:

```text
bed = empty
```

Multiple analyses may reference the same FASTQ paths.

---

# 30. Snakemake command

After finalization show the exact command required to run the pipeline.

For example:

```bash
snakemake \
    --snakefile /path/to/Snakefile \
    --configfile /path/to/run/config.yaml \
    --cores 16
```

Derive the actual command from the way the current pipeline is invoked.

Provide a copy-friendly presentation.

Do not run it automatically in the MVP.

---

# 31. Testing

Use tiny fixtures only.

Required tests include:

```text
SampleSheet parsing
FASTQ pairing
FASTQ pre-flight validation
upload-plan validation
filename safety
tusd hook handling
reference validation
primer validation
trimmed/untrimmed analyses
multiple analyses per FASTQ pair
analysis ID generation
tool option types
tool option range validation
unknown tool option rejection
raw command-line input rejection
run override resolution
analysis override resolution
samples.csv generation
config.yaml generation
Snakemake command generation
```

Also test the precedence:

```text
analysis override
    > run override
    > pipeline default
```

---

# 32. Explicitly out of scope

Do not implement:

```text
automatic Snakemake execution
pipeline progress monitoring
pipeline cancellation
results browser
authentication
user accounts
database
job queue
arbitrary shell commands
arbitrary CLI argument strings
automatic parsing of CLI help pages
React
Vue
Angular
```

These may be future additions.

---

# 33. Implementation order

Because development time is limited, implement in this order:

1. Inspect existing sample/config structure and tool invocations.
2. Define uploaded sample and analysis models.
3. Add reference/primer catalog.
4. Implement SampleSheet parsing.
5. Implement FASTQ filename matching.
6. Implement tool-option catalog and typed configuration models.
7. Add run/analysis configuration resolver to Snakemake.
8. Generate `samples.csv`.
9. Generate `config.yaml`.
10. Implement run-directory management.
11. Implement FastAPI basic application.
12. Implement SampleSheet upload.
13. Implement FASTQ upload plan and validation.
14. Integrate Uppy/tusd.
15. Implement analysis setup UI.
16. Implement run-wide advanced configuration UI.
17. Implement per-analysis advanced configuration UI.
18. Implement final validation.
19. Display Snakemake command.
20. Run the full relevant test suite.

Do not stop at architectural scaffolding.

Aim to have the complete end-to-end MVP functional before polishing secondary UI details.
