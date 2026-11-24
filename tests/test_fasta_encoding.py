"""FASTA guardados en Windows (latin-1, cp1252 o UTF-16) deben funcionar.

Regresión: un FASTA con caracteres no UTF-8 (p. ej. pegado desde Word)
hacía fallar el análisis con UnicodeDecodeError en Biopython.
"""
import pytest

from nanodelfy.core.refs import (
    build_edited_reference,
    decode_fasta_bytes,
    parse_first_fasta_record,
)

SEQ = "ACGT" * 25  # 100 nt


def _fasta_bytes(seq: str, header: str, encoding: str) -> bytes:
    text = f">{header}\n{seq}\n"
    return text.encode(encoding)


def test_decode_latin1_header():
    data = _fasta_bytes(SEQ, "referencia_ñ", "latin-1")
    assert "referencia_ñ" in decode_fasta_bytes(data)


def test_decode_utf16_with_bom():
    data = _fasta_bytes(SEQ, "referencia", "utf-16")
    assert decode_fasta_bytes(data).startswith(">referencia")


def test_decode_utf8_passthrough():
    data = _fasta_bytes(SEQ, "referencia", "utf-8")
    assert decode_fasta_bytes(data) == f">referencia\n{SEQ}\n"


def test_parse_first_record_from_latin1_file(tmp_path):
    fa = tmp_path / "wt_latin1.fa"
    fa.write_bytes(_fasta_bytes(SEQ, "WT_MmC7_ñ", "latin-1"))
    rec = parse_first_fasta_record(fa)
    assert str(rec.seq) == SEQ


def test_build_edited_reference_with_latin1_fasta(tmp_path):
    wt = tmp_path / "wt.fa"
    wt.write_bytes(_fasta_bytes(SEQ, "WT_áéíóú", "latin-1"))
    edit = tmp_path / "edit.fa"
    stats = build_edited_reference(wt, edit, deletion_start=11, deletion_end=20)
    assert stats["wt_length"] == 100
    assert stats["deletion_size"] == 10
    assert stats["edit_length"] == 90


def test_non_fasta_file_gives_clear_error(tmp_path):
    bad = tmp_path / "not_fasta.txt"
    bad.write_bytes(b"\xd2\x43\x47\x54 texto sin cabecera")
    with pytest.raises(ValueError, match="no parece un FASTA"):
        parse_first_fasta_record(bad)
