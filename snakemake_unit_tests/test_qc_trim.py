import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_qc_trim():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/qc_trim/data")
        expected_path = PurePosixPath(".tests/unit/qc_trim/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("/mnt/common/researchers/udo_gieraths/galaxy_port_data/results/1545554-HSV2_S82_not_trimmed_downsampled/qc/1545554-HSV2_S82_not_trimmed_downsampled_R1.fastq.gz /mnt/common/researchers/udo_gieraths/galaxy_port_data/results/1545554-HSV2_S82_not_trimmed_downsampled/qc/1545554-HSV2_S82_not_trimmed_downsampled_R2.fastq.gz /mnt/common/researchers/udo_gieraths/galaxy_port_data/results/1545554-HSV2_S82_not_trimmed_downsampled/qc/1545554-HSV2_S82_not_trimmed_downsampled.html /mnt/common/researchers/udo_gieraths/galaxy_port_data/results/1545554-HSV2_S82_not_trimmed_downsampled/qc/1545554-HSV2_S82_not_trimmed_downsampled.json", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "/mnt/common/researchers/udo_gieraths/galaxy_port_data/results/1545554-HSV2_S82_not_trimmed_downsampled/qc/1545554-HSV2_S82_not_trimmed_downsampled_R1.fastq.gz /mnt/common/researchers/udo_gieraths/galaxy_port_data/results/1545554-HSV2_S82_not_trimmed_downsampled/qc/1545554-HSV2_S82_not_trimmed_downsampled_R2.fastq.gz /mnt/common/researchers/udo_gieraths/galaxy_port_data/results/1545554-HSV2_S82_not_trimmed_downsampled/qc/1545554-HSV2_S82_not_trimmed_downsampled.html /mnt/common/researchers/udo_gieraths/galaxy_port_data/results/1545554-HSV2_S82_not_trimmed_downsampled/qc/1545554-HSV2_S82_not_trimmed_downsampled.json",
            "-f", 
            "-j1",
            "--keep-target-files",
            "--configfile",
            /mnt/common/researchers/udo_gieraths/code/CVG_Amplicon_Flow/my_config.yaml
    
            "--use-conda",
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
