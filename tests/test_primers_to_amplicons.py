from pathlib import Path

import pytest

from scripts.primers_to_amplicons import (
	read_primer_bed,
	build_amplicons,
	convert_primer_bed,
)


@pytest.fixture
def primer_bed(tmp_path: Path) -> Path:
	path = tmp_path / "primers.bed"
	path.write_text(
		"HSV_1_UL23\t0\t18\tUL23_0_LEFT\t1\t+\n"
		"HSV_1_UL23\t485\t506\tUL23_0_RIGHT\t1\t-\n"
		"HSV_1_UL23\t281\t301\tUL23_1_LEFT\t2\t+\n"
		"HSV_1_UL23\t840\t861\tUL23_1_RIGHT\t2\t-\n"
		"HSV_1_UL30\t0\t19\tUL30_0_LEFT\t1\t+\n"
		"HSV_1_UL30\t402\t426\tUL30_0_RIGHT\t1\t-\n"
	)
	return path


def test_build_amplicons_merges_primer_pairs() -> None:
	entries = read_primer_bed(
		Path("raw_data/HSV_1_UL23_UL30.bed")
	)

	amplicons = build_amplicons(entries)

	assert len(amplicons) == len(entries) // 2


def test_amplicon_spans_both_primers(primer_bed: Path) -> None:
	amplicons = build_amplicons(read_primer_bed(primer_bed))

	assert (amplicons[0].reference, amplicons[0].start, amplicons[0].end,
			amplicons[0].name) == ("HSV_1_UL23", 0, 506, "UL23_0")
	assert (amplicons[1].reference, amplicons[1].start, amplicons[1].end,
			amplicons[1].name) == ("HSV_1_UL23", 281, 861, "UL23_1")
	assert (amplicons[2].reference, amplicons[2].start, amplicons[2].end,
			amplicons[2].name) == ("HSV_1_UL30", 0, 426, "UL30_0")


def test_amplicons_have_score_one_and_plus_strand(primer_bed: Path) -> None:
	amplicons = build_amplicons(read_primer_bed(primer_bed))

	for amplicon in amplicons:
		assert amplicon.score == 1
		assert amplicon.strand == "+"


def test_convert_primer_bed_writes_bed_file(primer_bed: Path, tmp_path: Path) -> None:
	output = tmp_path / "amplicons.bed"

	convert_primer_bed(primer_bed, output)

	lines = output.read_text().splitlines()
	assert lines == [
		"HSV_1_UL23\t0\t506\tUL23_0\t1\t+",
		"HSV_1_UL23\t281\t861\tUL23_1\t1\t+",
		"HSV_1_UL30\t0\t426\tUL30_0\t1\t+",
	]

def test_odd_number_of_entries_raises(tmp_path: Path) -> None:
	path = tmp_path / "odd.bed"
	path.write_text(
		"ref\t0\t18\tA_0_LEFT\t1\t+\n"
		"ref\t485\t506\tA_0_RIGHT\t1\t-\n"
		"ref\t900\t910\tA_1_LEFT\t1\t+\n"
	)

	with pytest.raises(ValueError, match="even number of entries"):
		build_amplicons(read_primer_bed(path))