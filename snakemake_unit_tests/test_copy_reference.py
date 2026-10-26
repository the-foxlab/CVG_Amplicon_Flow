import os
import csv
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_copy_reference():

    with TemporaryDirectory() as tmpdir:
        test_dir = Path(__file__).resolve().parent
        workdir = Path(tmpdir) / "workdir"
        data_path = test_dir / "copy_reference" / "data"
        expected_path = test_dir / "copy_reference" / "expected"
        target = (
            workdir
            / "results"
            / "references"
            / "83985b3dd77c1fe33ef74fc2b3847ae38467da9dd323cf3d1ab12d07957cc021"
            / "83985b3dd77c1fe33ef74fc2b3847ae38467da9dd323cf3d1ab12d07957cc021.unmasked.fasta"
        )

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        source_sheet = test_dir / "sample_sheet.csv"
        sample_sheet = workdir / "sample_sheet.csv"
        with source_sheet.open(newline="") as input_handle:
            reader = csv.DictReader(input_handle)
            fieldnames = reader.fieldnames
            rows = list(reader)

        if fieldnames is None:
            raise ValueError(f"Sample sheet has no header: {source_sheet}")

        input_columns = ("reference_fasta", "bed", "r1", "r2")
        for row in rows:
            for column in input_columns:
                value = row[column].strip()
                if value:
                    row[column] = str((test_dir / value).resolve())

        with sample_sheet.open("w", newline="") as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

        # dbg
        print(target, file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            str(target),
            "-f", 
            "-j1",
            "--keep-target-files",
            "--allowed-rules",
            "copy_reference",
            "--configfile",
            f"{test_dir}/test_config.yaml",
            "--use-conda",
            "--config",
            f"sample_sheet={sample_sheet}",
            f"res_dir={workdir / 'results'}",
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
