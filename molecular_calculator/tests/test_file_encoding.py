"""Tests for robust CSV encoding handling on file upload.

Regression coverage for the bug where uploading a CSV exported from Excel or
a Windows-based tool raised:

    Error reading file: 'utf-8' codec can't decode byte 0xb5 ... invalid start byte

The 0xb5 byte is the micro sign (``µ``) in Windows-1252/Latin-1, commonly
present in potency unit headers such as ``Potency µM``. ``pd.read_csv``
defaults to UTF-8, so such files failed to load.
"""

from io import BytesIO

import pandas as pd
import pytest

from molecular_calculator.ui.components.file_uploader import (
    read_csv_robust,
    read_uploaded_file,
)


class _FakeUpload(BytesIO):
    """Minimal stand-in for a Streamlit UploadedFile (BytesIO + ``name``)."""

    def __init__(self, content: bytes, name: str):
        super().__init__(content)
        self.name = name


# A CSV whose header contains the micro sign (0xb5 in cp1252) -> not valid UTF-8.
_CP1252_CSV = "SMILES,Potency µM\nCCO,5\nc1ccccc1,10\n".encode("cp1252")
_UTF8_CSV = "SMILES,Potency µM\nCCO,5\nc1ccccc1,10\n".encode("utf-8")
_UTF8_BOM_CSV = "SMILES,Name\nCCO,Ethanol\n".encode("utf-8-sig")


def test_plain_read_csv_reproduces_bug():
    """Sanity check: the default UTF-8 read is what fails on cp1252 input."""
    with pytest.raises(UnicodeDecodeError):
        pd.read_csv(BytesIO(_CP1252_CSV))


def test_read_csv_robust_handles_cp1252():
    df = read_csv_robust(BytesIO(_CP1252_CSV))
    assert list(df.columns) == ["SMILES", "Potency µM"]
    assert len(df) == 2


def test_read_csv_robust_handles_plain_utf8():
    df = read_csv_robust(BytesIO(_UTF8_CSV))
    assert list(df.columns) == ["SMILES", "Potency µM"]


def test_read_csv_robust_strips_utf8_bom():
    df = read_csv_robust(BytesIO(_UTF8_BOM_CSV))
    # Without BOM stripping the first column name would be '﻿SMILES'.
    assert list(df.columns) == ["SMILES", "Name"]


def test_read_uploaded_file_reads_cp1252_csv():
    upload = _FakeUpload(_CP1252_CSV, "compounds.csv")
    df = read_uploaded_file(upload, show_errors=False)
    assert df is not None
    assert "Potency µM" in df.columns


def test_bom_export_round_trips_and_strips_bom():
    """Our CSV exports use utf-8-sig (BOM for Excel); re-import must be clean."""
    df = pd.DataFrame({"SMILES": ["CCO"], "MIC min (µM)": [5]})
    exported = df.to_csv(index=False).encode("utf-8-sig")

    # The BOM must be present so Excel auto-detects UTF-8.
    assert exported[:3] == b"\xef\xbb\xbf"

    # Re-importing our own export strips the BOM (no '﻿SMILES' column).
    back = read_csv_robust(BytesIO(exported))
    assert list(back.columns) == ["SMILES", "MIC min (µM)"]


def test_utf16_csv_is_read_correctly_not_as_garbage():
    raw = "SMILES,val\nCCO,1\nc1ccccc1,2\n".encode("utf-16")  # has BOM 0xFF 0xFE
    df = read_csv_robust(BytesIO(raw))
    assert list(df.columns) == ["SMILES", "val"]
    assert df["SMILES"].tolist() == ["CCO", "c1ccccc1"]


def test_nul_bytes_in_decoded_header_raise_clear_error():
    # UTF-16 content WITHOUT a usable BOM, forced through latin-1, contains NULs.
    raw = "A,B\n1,2\n".encode("utf-16-le")  # no BOM
    with pytest.raises(Exception) as exc:
        read_csv_robust(BytesIO(raw))
    assert "encod" in str(exc.value).lower()
