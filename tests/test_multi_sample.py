"""Varios FASTQ en la misma petición: una muestra por archivo.

El endpoint /api/run acepta múltiples archivos bajo el campo "fastq";
cada uno se analiza con las mismas condiciones y el trabajo termina con
result["samples"] = [{name, result}, ...], donde el nombre sale del
archivo. Las descargas llevan el prefijo s<índice>_.
"""
import random
import time
from pathlib import Path

import pytest

pytest.importorskip("httpx", reason="TestClient requiere httpx")

from fastapi.testclient import TestClient

from nanodelfy.server import app


def _datos(tmp_path: Path, prefijo: str):
    """Amplicón WT pequeño con una deleción de 200 nt y lecturas editadas."""
    rng = random.Random(21)
    wt = list("".join(rng.choice("ACGT") for _ in range(700)))
    wt = "".join(wt)
    (tmp_path / "wt.fa").write_text(">amp\n" + wt + "\n")
    plantilla = wt[330:353]
    ins = plantilla[:10] + ("A" if plantilla[10] != "A" else "C") + plantilla[11:]
    reads = [(f"{prefijo}_edit_{i}", wt[:199] + ins + wt[400:])
             for i in range(8)]
    reads += [(f"{prefijo}_wt_{i}", wt) for i in range(12)]
    fq = tmp_path / f"{prefijo}.fastq"
    with open(fq, "w") as f:
        for rid, s in reads:
            f.write(f"@{rid}\n{s}\n+\n{'I' * len(s)}\n")
    return fq


def test_varios_fastq_en_un_trabajo(tmp_path):
    client = TestClient(app)
    fq1 = _datos(tmp_path, "Phire_MEF5")
    fq2 = _datos(tmp_path, "Phire_MEF8")

    with open(fq1, "rb") as f1, open(fq2, "rb") as f2, \
            open(tmp_path / "wt.fa", "rb") as fw:
        r = client.post("/api/run", files=[
            ("fastq", ("Phire_MEF5_270526.fastq", f1, "application/octet-stream")),
            ("fastq", ("Phire_MEF8_270526-2.fastq", f2, "application/octet-stream")),
            ("wt_fasta", ("wt.fa", fw, "text/plain")),
        ], data={
            "deletion_start": "200", "deletion_end": "400",
            "min_quality": "10", "min_length": "400", "max_length": "900",
            "breakpoint_tolerance": "20", "min_sv_deletion": "100",
            "nhej_search_window": "200", "wt_min": "650", "wt_max": "760",
        })
    assert r.status_code == 200, r.text
    job_id = r.json()["job_id"]

    # esperar a que termine (el worker corre en un hilo aparte)
    deadline = time.time() + 180
    job = None
    while time.time() < deadline:
        job = client.get(f"/api/job/{job_id}").json()
        if job["status"] in ("done", "error"):
            break
        time.sleep(1)
    assert job["status"] == "done", "\n".join(job["log"][-20:])

    muestras = job["result"]["samples"]
    nombres = [m["name"] for m in muestras]
    assert nombres == ["Phire_MEF5_270526", "Phire_MEF8_270526-2"], nombres
    for m in muestras:
        res = m["result"]
        assert res["classification"]["summary"], "sin resumen de clasificación"
        assert res["deletion_view"]["deleciones"], "sin deleciones recurrentes"

    # descargas por muestra (prefijo s<índice>_) y zip global
    r = client.get(f"/api/download/{job_id}/s0_excel")
    assert r.status_code == 200 and r.headers["content-type"].startswith(
        "application/vnd.openxmlformats"), r.status_code
    r = client.get(f"/api/download/{job_id}/s1_summary_tsv")
    assert r.status_code == 200, r.status_code
    r = client.get(f"/api/download/{job_id}/all")
    assert r.status_code == 200 and r.headers["content-type"] == "application/zip"
    # sin prefijo no es válido en el formato nuevo
    r = client.get(f"/api/download/{job_id}/excel")
    assert r.status_code == 404


def test_un_solo_fastq_sigue_funcionando(tmp_path):
    """Un único archivo: mismo formato nuevo, etiqueta = nombre del archivo."""
    client = TestClient(app)
    fq = _datos(tmp_path, "muestra_unica")

    with open(fq, "rb") as f, open(tmp_path / "wt.fa", "rb") as fw:
        r = client.post("/api/run", files=[
            ("fastq", ("muestra_unica.fastq", f, "application/octet-stream")),
            ("wt_fasta", ("wt.fa", fw, "text/plain")),
        ], data={
            "deletion_start": "200", "deletion_end": "400",
            "min_quality": "10", "min_length": "400", "max_length": "900",
            "breakpoint_tolerance": "20", "min_sv_deletion": "100",
            "nhej_search_window": "200", "wt_min": "650", "wt_max": "760",
        })
    assert r.status_code == 200, r.text
    job_id = r.json()["job_id"]

    deadline = time.time() + 180
    job = None
    while time.time() < deadline:
        job = client.get(f"/api/job/{job_id}").json()
        if job["status"] in ("done", "error"):
            break
        time.sleep(1)
    assert job["status"] == "done", "\n".join(job["log"][-20:])

    # etiqueta por defecto = nombre de la muestra (del archivo)
    assert job["sample_name"] == "muestra_unica", job["sample_name"]
    muestras = job["result"]["samples"]
    assert len(muestras) == 1 and muestras[0]["name"] == "muestra_unica"
    r = client.get(f"/api/download/{job_id}/s0_summary_tsv")
    assert r.status_code == 200
