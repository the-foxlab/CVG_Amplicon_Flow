#!/usr/bin/env bash
set -euo pipefail

: "${GALAXY_PORT_TUSD_UPLOAD_ROOT:=/tmp/galaxy_port_tusd}"
: "${GALAXY_PORT_TUSD_HOOK_SECRET:=dev-hook-secret}"
: "${GALAXY_PORT_FASTAPI_BASE:=http://127.0.0.1:8000}"

mkdir -p "${GALAXY_PORT_TUSD_UPLOAD_ROOT}"

exec tusd \
	-host 127.0.0.1 \
	-port 1080 \
	-base-path /files/ \
	-upload-dir "${GALAXY_PORT_TUSD_UPLOAD_ROOT}" \
	-hooks-http "${GALAXY_PORT_FASTAPI_BASE}/internal/tusd/hooks" \
	-hooks-enabled-events pre-create,post-finish \
	-hooks-http-forward-headers X-Galaxy-Port-Hook-Secret
