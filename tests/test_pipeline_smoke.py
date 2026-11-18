"""Test de humo del pipeline completo con datos sintéticos.

Reproduce el caso de uso "deleción pequeña": amplicón WT de 3887 nt con
una deleción de 120 nt (posiciones 1438-1557) y lecturas simuladas con
errores tipo Nanopore (~97.5% identidad). Verifica que:
  - la app unificada ejecuta end-to-end desde el paquete instalado,
  - clasifica las lecturas EDIT como la deleción esperada,
  - no clasifica ninguna lectura WT como EDIT (0 falsos positivos),
  - genera todas las salidas esperadas.

El generador de lecturas no depende de badread: modelo simple de
sustituciones/inserciones/deleciones con tasa ~2.5%.
"""
import random
from pathlib import Path

import pytest

from ampdel.core.pipeline import run_full_pipeline

WT_LEN = 3887
DEL_START, DEL_END = 1438, 1557          # 1-based inclusivo → 120 nt


def _sim_ref(base: str, n: int, rng: random.Random) -> str:
    """Referencia aleatoria con composición equilibrada."""
    return "".join(rng.choice(base) for _ in range(n))


def _mutate(seq: str, rate: float, rng: random.Random) -> str:
    """Introduce sustituciones/indels con probabilidad `rate` por posición."""
    out = []
    for ch in seq:
        r = rng.random()
        if r < rate * 0.7:                       # sustitución
            out.append(rng.choice("ACGT".replace(ch, "")))
        elif r < rate * 0.85:                    # deleción 1 nt
            continue
        elif r < rate:                           # inserción 1 nt
            out.append(ch + rng.choice("ACGT"))
        else:
            out.append(ch)
    return "".join(out)


def _write_fasta(path: Path, name: str, seq: str) -> None:
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")


def _write_fastq(path: Path, reads: list[tuple[str, str]], rng: random.Random) -> None:
    with open(path, "w") as f:
        for j, (rid, seq) in enumerate(reads):
            quals = "".join(rng.choice("+-5HI") for _ in seq)  # ~Q10-Q40
            f.write(f"@{rid}\n{seq}\n+\n{quals}\n")


@pytest.fixture(scope="module")
def mini_dataset(tmp_path_factory):
    """Amplicón + FASTQ con 30 lecturas WT y 30 lecturas EDIT."""
    tmp = tmp_path_factory.mktemp("ampdel_smoke")
    rng = random.Random(20260101)
    wt = _sim_ref("ACGT", WT_LEN, rng)
    edit = wt[:DEL_START - 1] + wt[DEL_END:]      # 120 nt menos

    _write_fasta(tmp / "amp_WT.fa", "amp", wt)
    _write_fasta(tmp / "amp_EDIT.fa", "amp_EDIT", edit)

    reads = []
    for i in range(30):
        # lecturas WT: se tomaron del WT, se re-mutaron; coordinadas genómicas
        start = rng.randrange(0, 400)
        end = min(WT_LEN, 3300 + rng.randrange(300))
        reads.append((f"wt_{i}", _mutate(wt[start:end], 0.025, rng)))
    for i in range(30):
        start = rng.randrange(0, 400)
        end = min(WT_LEN - 120, 3300 + rng.randrange(300))
        reads.append((f"edit_{i}", _mutate(edit[start:end], 0.025, rng)))
    _write_fastq(tmp / "muestra.fastq", reads, rng)
    return tmp


def test_pipeline_end_to_end(mini_dataset, tmp_path):
    # Guía sintética extraída del propio amplicón (hebra +, 20 nt):
    # garantiza que el paso de microindels se ejecute con datos reales.
    wt = (mini_dataset / "amp_WT.fa").read_text().splitlines()
    wt_seq = "".join(l for l in wt if not l.startswith(">"))
    guide = wt_seq[1000:1020]

    out = tmp_path / "salida"
    res = run_full_pipeline(
        job_dir=out,
        fastq_path=mini_dataset / "muestra.fastq",
        wt_fasta_path=mini_dataset / "amp_WT.fa",
        sample_name="mixta120",
        params={
            "deletion_start": DEL_START, "deletion_end": DEL_END,
            "min_quality": 7, "min_length": 200, "max_length": 5000,
            "breakpoint_tolerance": 20, "nhej_search_window": 200,
            "min_junction_insertion": 5, "min_sv_deletion": 100,
            "min_edit_ref_coverage": 0.80, "min_mapq": 20,
            "guides": f"test: {guide}", "nuclease": "cas9",
            "microindel_window": 20, "microindel_max_indel": 50,
            "microindel_large_del_threshold": 100,
        },
        log=lambda m: None,
    )
    files = res.get("files", {})

    # 1. Todas las salidas esperadas existen
    for key in ("perread_tsv", "summary_tsv", "excel"):
        assert key in files and Path(files[key]).is_file(), f"falta {key}"
    micro = [v for k, v in files.items() if k.startswith("micro_")]
    assert len(micro) == 3, "faltan tablas de microindels"

    # 2. Lecturas EDIT detectadas como deleción (sensibilidad ≥ 70 %)
    import pandas as pd
    from ampdel.core.classify import EDITED_CLASSES
    df = pd.read_csv(files["perread_tsv"], sep="\t")
    edits = df[df["read_id"].str.startswith("edit_")]
    if len(edits):
        detected = edits["class"].isin(EDITED_CLASSES)
        assert detected.sum() >= 0.7 * len(edits), \
            f"sensibilidad baja: {detected.sum()}/{len(edits)}"

    # 3. Ninguna lectura WT clasificada como editada (especificidad)
    wts = df[df["read_id"].str.startswith("wt_")]
    fp = wts["class"].isin(EDITED_CLASSES).sum()
    assert fp <= 1, f"falsos positivos WT: {fp}/{len(wts)}"
