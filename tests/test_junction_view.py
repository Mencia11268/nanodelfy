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


def test_variantes_de_borde_de_insercion_se_agrupan(tmp_path):
    """Inserciones que difieren en 1 nt del borde son la misma variante.

    El ruido de homopolímeros de Nanopore (p. ej. 3A que pasan a 2) hace que
    el alineador deslice el anclaje de una inserción y muestre bases de
    compensación en el borde 5' que no son inserciones reales: sin la
    agrupación por similitud, la misma variante aparece como varios grupos
    de pocas lecturas. Aquí las tres formas (exacta, con una base extra en
    5' y sin la última base) deben quedar en UN grupo cuya inserción
    mostrada es la forma más frecuente.
    """
    rng = random.Random(42)
    wt = "".join(rng.choice("ACGT") for _ in range(WT_LEN))
    _write_fasta(tmp_path / "wt.fa", wt)

    plantilla = wt[1300:1320]
    ins = plantilla[:10] + ("A" if plantilla[10] != "A" else "C") + plantilla[11:]
    formas = [("exacta", ins, 10), ("con_A_5", "A" + ins, 5),
              ("sin_1_3", ins[:-1], 5)]

    reads = []
    for nombre, forma, n in formas:
        reads += [(f"ins_{nombre}_{i}", wt[:DS - 1] + forma + wt[DE:])
                  for i in range(n)]
    reads += [(f"precise_{i}", wt[:DS - 1] + wt[DE:]) for i in range(12)]
    reads += [(f"wt_{i}", wt) for i in range(20)]
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

    delecs = res["deletion_view"]["deleciones"]
    con_ins = [d for d in delecs if d["ins_union"]]
    assert len(con_ins) == 1, \
        f"las tres formas deben agruparse en una, hay {len(con_ins)}"
    assert con_ins[0]["n"] == 20, \
        f"el grupo debe reunir las 20 lecturas: {con_ins[0]['n']}"
    assert con_ins[0]["ins_union"] == ins, \
        f"debe mostrarse la forma más frecuente: {con_ins[0]['ins_union']!r}"
