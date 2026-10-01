import logging
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from scripts.mask_refs import (
	Primer,
	mask_fasta_beginnings,
	prepare_fastas,
	read_primers,
)


@pytest.fixture
def primer_bed(tmp_path: Path) -> Path:
	path = tmp_path / "primers.bed"
	path.write_text(
		"# reference\tstart\tend\n"
		"ref1\t12\t15\n"
		"ref1\t4\t6\n"
		"ref2\t0\t2\n"
		"\n"
	)
	return path


@pytest.fixture
def fasta_file(tmp_path: Path) -> Path:
	path = tmp_path / "references.fasta"
	path.write_text(
		">ref1 first reference\n"
		"ACGTACGTACGTACGTACGT\n"
		">ref3 without primers\n"
		"TTGCA\n"
	)
	return path


def test_read_primers_skips_comments_and_blank_lines(primer_bed: Path) -> None:
	assert read_primers(primer_bed) == {
		"ref1": [Primer("ref1", 12, 15), Primer("ref1", 4, 6)],
		"ref2": [Primer("ref2", 0, 2)],
	}


@pytest.mark.parametrize(
	("bed_content", "message"),
	[
		("ref1\t2\n", "expected at least 3 columns"),
		("ref1\tstart\t5\n", "Invalid BED coordinates on line 1"),
		("ref1\t-1\t5\n", "Invalid BED interval on line 1"),
		("ref1\t5\t5\n", "Invalid BED interval on line 1"),
	],
)
def test_read_primers_rejects_invalid_rows(
	tmp_path: Path,
	bed_content: str,
	message: str,
) -> None:
	path = tmp_path / "invalid.bed"
	path.write_text(bed_content)

	with pytest.raises(ValueError, match=message):
		read_primers(path)


def test_mask_fasta_beginnings_masks_outer_ends() -> None:
	record = SeqRecord(Seq("ACGTACGTACGTACGTACGT"), id="ref1")
	primers = [Primer("ref1", 12, 15), Primer("ref1", 4, 6)]

	result = mask_fasta_beginnings(record, primers)

	assert result is record
	assert str(result.seq) == "NNNNNNGTACGTNNNNNNNN"


def test_mask_fasta_beginnings_without_primers_leaves_record_unchanged() -> None:
	record = SeqRecord(Seq("ACGT"), id="ref1")

	result = mask_fasta_beginnings(record, [])

	assert result is record
	assert str(result.seq) == "ACGT"


@pytest.mark.parametrize(
	"primers",
	[
		[Primer("ref1", 0, 21)],
		[Primer("ref1", 0, 1), Primer("ref1", 21, 22)],
	],
)
def test_mask_fasta_beginnings_rejects_intervals_past_sequence(
	primers: list[Primer],
) -> None:
	record = SeqRecord(Seq("A" * 20), id="ref1")

	with pytest.raises(ValueError, match="exceeds sequence length"):
		mask_fasta_beginnings(record, primers)


def test_prepare_fastas_masks_matching_records_and_warns_for_unmatched(
	fasta_file: Path,
	primer_bed: Path,
	caplog: pytest.LogCaptureFixture,
) -> None:
	with caplog.at_level(logging.WARNING):
		records = list(prepare_fastas(fasta_file, primer_bed))

	assert [record.id for record in records] == ["ref1", "ref3"]
	assert str(records[0].seq) == "NNNNNNGTACGTNNNNNNNN"
	assert str(records[1].seq) == "TTGCA"
	assert "No primers found for reference 'ref3'" in caplog.text