"""Module entrypoint for local development."""

from __future__ import annotations

import uvicorn


def main() -> None:
	"""Run the FastAPI application with uvicorn."""

	uvicorn.run("galaxy_port.webapp:app", host="127.0.0.1", port=8000, reload=True)


if __name__ == "__main__":
	main()
