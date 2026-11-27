"""Vista de uniones: deleciones con inserciones NHEJ y deleciones de flanco.

Escenario: amplicón WT de 3000 nt con la deleción esperada [1000, 1500)
(1-based 1000-1499... se pasa como deletion_start=1000, deletion_end=1500).
Tres poblaciones de lecturas editadas:

  - precisa:      wt[:999] + wt[1500:]
  - con inserción: wt[:999] + "GATTACAGCCTA" + wt[1500:]   (NHEJ, 12 nt)
  - con resección 5': wt[:980] + wt[985:999] + wt[1500:]  (deleción de 5 nt
    en el flanco 5' SEPARADA de la unión por bases coincidentes: las 5 bases
    del flanco van tachadas y la barra no cambia)

El gráfico de la web (deletion_view) debe separar las tres uniones y
representar correctamente inserción de unión y deleción de flanco.
"""
import random
from pathlib import Path

from nanodelfy.core.pipeline import run_full_pipeline

WT_LEN = 3000
DS, DE = 1000, 1500          # 1-based: ds = 999, de = 1500 (0-based semiabierto)


def _write_fasta(path: Path, seq: str) -> None:
    with open(path, "w") as f:
        f.write(">amp\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")


def test_vista_de_union_con_insercion_y_reseccion(tmp_path):
    rng = random.Random(42)
    wt = "".join(rng.choice("ACGT") for _ in range(WT_LEN))
    _write_fasta(tmp_path / "wt.fa", wt)

    precise = wt[:DS - 1] + wt[DE:]
    with_ins = wt[:DS - 1] + "GATTACAGCCTA" + wt[DE:]
    with_res = wt[:DS - 20] + wt[DS - 15:DS - 1] + wt[DE:]  # del 5 nt en flanco 5'

    reads = ([(f"precise_{i}", precise) for i in range(20)]
             + [(f"ins_{i}", with_ins) for i in range(10)]
             + [(f"res_{i}", with_res) for i in range(8)]
             + [(f"wt_{i}", wt) for i in range(20)])
    rng.shuffle(reads)
    fq = tmp_path / "muestra.fastq"
    with open(fq, "w") as f:
        for rid, s in reads:
            f.write(f"@{rid}\n{s}\n+\n{'I' * len(s)}\n")

    res = run_full_pipeline(
        tmp_path / "job", fq, tmp_path / "wt.fa", "mix",
        params={
            "deletion_start": DS, "deletion_end": DE,
            "min_quality": 10, "min_length": 1500, "max_length": 3500,
            "breakpoint_tolerance": 20, "min_sv_deletion": 100,
            "nhej_search_window": 200, "wt_min": 2900, "wt_max": 3100,
        },
        log=lambda m: None,
    )

    dv = res["deletion_view"]
    assert dv, "deletion_view debe existir"
    delecs = dv["deleciones"]
    assert delecs, "debe haber deleciones recurrentes"
    # tres uniones distintas, ordenadas por frecuencia
    assert len(delecs) >= 3, f"esperaba 3 uniones, hay {len(delecs)}"
    top = delecs[0]
    assert top["n"] == 20, f"la unión precisa debe agrupar 20: {top['n']}"

    # unión con inserción NHEJ: la inserción de 12 nt se reparte entre
    # inserción de unión (azul) e inserción anclada en el flanco 5'
    # (también azul): juntas deben recuperar la mayor parte de GATTACAGCCTA
    con_ins = [d for d in delecs
               if d["ins_union"]
               or any(s["t"] == "i" for s in d["flanco5"])]
    assert con_ins, "debe detectarse la unión con inserción"
    assert con_ins[0]["n"] == 10
    azules = con_ins[0]["ins_union"] + "".join(
        s["seq"] for s in con_ins[0]["flanco5"] if s["t"] == "i")
    assert len(azules) >= 9, \
        f"la inserción NHEJ debe verse en azul, hay {azules!r}"

    # unión con resección 5': 5 bases tachadas en el flanco 5'
    con_res = [d for d in delecs
               if any(s["t"] == "d" for s in d["flanco5"])]
    assert con_res, "debe detectarse la unión con deleción en el flanco 5'"
    r = con_res[0]
    assert r["n"] == 8
    d_segs = [s for s in r["flanco5"] if s["t"] == "d"]
    assert sum(len(s["seq"]) for s in d_segs) == 5, \
        "la resección de 5 nt debe quedar tachada en el flanco 5'"
    assert wt[DS - 20:DS - 15] == "".join(s["seq"] for s in d_segs), \
        "las bases tachadas deben ser las de la resección"
    # los segmentos m+d del flanco cubren los 30 nt de referencia
    cubierto = sum(len(s["seq"]) for s in r["flanco5"] if s["t"] != "i")
    assert cubierto == 30


def _guia_cas9_mas_arriba(wt: str, desde: int) -> tuple:
    """Primer protospacer '+' con PAM NGG cuyo inicio sea >= desde.

    Devuelve (inicio, spacer): el corte romo cae en inicio + 17.
    """
    for a in range(desde, len(wt) - 23):
        pam = wt[a + 20:a + 23]
        if pam[1:] == "GG":
            return a, wt[a:a + 20]
    raise AssertionError("no hay guía con PAM NGG en el tramo")


def test_deleccion_teorica_desde_guias(tmp_path):
    """Dos guías cas9 → deletion_view.teorica con el tamaño esperado.

    El corte romo de cas9 cae 3 nt aguas arriba del PAM (entre las bases
    17|18 del protospacer de 20 nt), así que la deleción teórica entre
    dos guías es la distancia entre sus cortes.
    """
    rng = random.Random(42)
    wt = "".join(rng.choice("ACGT") for _ in range(WT_LEN))
    _write_fasta(tmp_path / "wt.fa", wt)

    a1, sp1 = _guia_cas9_mas_arriba(wt, 950)
    a2, sp2 = _guia_cas9_mas_arriba(wt, 1400)
    corte1, corte2 = a1 + 17, a2 + 17
    assert corte1 < corte2

    precise = wt[:DS - 1] + wt[DE:]
    reads = ([(f"precise_{i}", precise) for i in range(20)]
             + [(f"wt_{i}", wt) for i in range(20)])
    rng.shuffle(reads)
    fq = tmp_path / "muestra2.fastq"
    with open(fq, "w") as f:
        for rid, s in reads:
            f.write(f"@{rid}\n{s}\n+\n{'I' * len(s)}\n")

    res = run_full_pipeline(
        tmp_path / "job2", fq, tmp_path / "wt.fa", "mix2",
        params={
            "deletion_start": DS, "deletion_end": DE,
            "min_quality": 10, "min_length": 1500, "max_length": 3500,
            "breakpoint_tolerance": 20, "min_sv_deletion": 100,
            "nhej_search_window": 200, "wt_min": 2900, "wt_max": 3100,
            "guides": f"guia_A: {sp1}\nguia_B: {sp2}",
            "nuclease": "cas9", "guide_max_mismatch": 0,
        },
        log=lambda m: None,
    )

    t = res["deletion_view"]["teorica"]
    assert t, "debe calcularse la deleción teórica"
    assert t["guia5"] == "guia_A" and t["guia3"] == "guia_B"
    assert t["nuclease"] == "cas9"
    assert t["start"] == corte1 and t["end"] == corte2
    assert t["size"] == corte2 - corte1
    # el resto de la vista no se ve afectado
    assert res["deletion_view"]["deleciones"][0]["n"] == 20
