# Local Development Stack

This project uses FastAPI as control plane and tusd as upload data plane.

## 1. Install Python dependencies

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements/requirements_gui_mvp.txt
```

## 2. Install pinned frontend dependencies

```bash
npm install
```

This installs local copies for:
- `@uppy/core`
- `@uppy/dashboard`
- `@uppy/tus`
- `htmx.org`

The application serves these from `node_modules` via `/vendor/...`.

## 3. Start FastAPI (localhost only)

```bash
export GALAXY_PORT_STORAGE_ROOT=/tmp/galaxy_port_runs
export GALAXY_PORT_TUSD_UPLOAD_ROOT=/tmp/galaxy_port_tusd
export GALAXY_PORT_TUSD_HOOK_SECRET=dev-hook-secret
.venv/bin/python -m galaxy_port
```

## 4. Start tusd (localhost only)

```bash
chmod +x deploy/dev/run_tusd.sh
GALAXY_PORT_TUSD_HOOK_SECRET=dev-hook-secret ./deploy/dev/run_tusd.sh
```

Notes:
- tusd listens on `127.0.0.1:1080`
- tus upload base path is `/files/`
- tus hooks call FastAPI at `/internal/tusd/hooks`
- FastAPI does not proxy FASTQ bytes

## 5. Start nginx as browser-facing reverse proxy

```bash
nginx -c "$PWD/deploy/dev/nginx.conf"
```

Open:
- `http://127.0.0.1:8088/`

Routing:
- `/` -> FastAPI (`127.0.0.1:8000`)
- `/files/` -> tusd (`127.0.0.1:1080/files/`)

The hook endpoint remains on FastAPI localhost and is not separately exposed.

## 6. Stop nginx

```bash
nginx -s stop
```

## 7. Typical local flow

1. Create staging run.
2. Upload `SampleSheet.csv`.
3. Select FASTQs in Uppy dashboard.
4. Click **Validate FASTQ Selection**.
5. Click **Start Upload**.
6. Add one or more analyses.
7. Optionally set run-level and per-analysis overrides.
8. Click **Prepare Run**.
9. Copy Snakemake command from prepared view.
