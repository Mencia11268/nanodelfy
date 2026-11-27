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


def test_coordenadas_de_deleccion_derivadas_de_guias(tmp_path):
    """Sin deletion_start/end: se derivan de los cortes de dos guías.

    Las lecturas precisas deben clasificarse DEL_EXACT con la deleción
    derivada de las guías, y `esperada` debe coincidir con `teorica`.
    """
    rng = random.Random(42)
    wt = "".join(rng.choice("ACGT") for _ in range(WT_LEN))
    _write_fasta(tmp_path / "wt.fa", wt)

    a1, sp1 = _guia_cas9_mas_arriba(wt, 950)
    a2, sp2 = _guia_cas9_mas_arriba(wt, 1400)
    corte1, corte2 = a1 + 17, a2 + 17
    # deleción que dejarían los dos cortes: 1-based [corte1+1, corte2]
    edited = wt[:corte1] + wt[corte2:]

    reads = ([(f"precise_{i}", edited) for i in range(20)]
             + [(f"wt_{i}", wt) for i in range(20)])
    rng.shuffle(reads)
    fq = tmp_path / "muestra3.fastq"
    with open(fq, "w") as f:
        for rid, s in reads:
            f.write(f"@{rid}\n{s}\n+\n{'I' * len(s)}\n")

    res = run_full_pipeline(
        tmp_path / "job3", fq, tmp_path / "wt.fa", "mix3",
        params={
            "deletion_start": None, "deletion_end": None,
            "min_quality": 10, "min_length": 1500, "max_length": 3500,
            "breakpoint_tolerance": 20, "min_sv_deletion": 100,
            "nhej_search_window": 200,
            "guides": f"guia_A: {sp1}\nguia_B: {sp2}",
            "nuclease": "cas9", "guide_max_mismatch": 0,
        },
        log=lambda m: None,
    )

    dv = res["deletion_view"]
    t = dv["teorica"]
    assert t and t["size"] == corte2 - corte1
    # las coordenadas derivadas pasan a ser la deleción esperada
    assert dv["esperada"]["start"] == corte1
    assert dv["esperada"]["end"] == corte2
    # las lecturas precisas (según los cortes) son DEL_EXACT
    resumen = {r["class"]: r["reads"]
               for r in res["classification"]["summary"]}
    assert resumen.get("DEL_EXACT") == 20, resumen


def test_recorte_de_bordes_templados_unitario():
    from nanodelfy.core.classify import _recorte_bordes_templados
    rng = random.Random(99)
    region = "".join(rng.choice("ACGT") for _ in range(200))
    core = region[50:73]                     # copia exacta: 23 nt
    core1 = core[:10] + ("A" if core[10] != "A" else "C") + core[11:]

    # copia limpia y copia con 1 discrepancia: intactas
    assert _recorte_bordes_templados(core, region) == core
    assert _recorte_bordes_templados(core1, region) == core1
    # bases extra en el borde 5' SIN pérdida en el flanco: son una
    # inserción real y se conservan (aunque la copia sola alinee mejor)
    assert _recorte_bordes_templados("G" + core, region) == "G" + core
    assert _recorte_bordes_templados("GC" + core, region) == "GC" + core
    # bases extra en el borde 3': igual
    assert _recorte_bordes_templados(core + "G", region) == core + "G"
    # bases extra CON pérdida equivalente en un homopolímero del flanco:
    # son compensación y se recortan
    assert _recorte_bordes_templados("G" + core, region, perdidas5=1) == core
    assert _recorte_bordes_templados("GC" + core, region, perdidas5=2) == core
    assert _recorte_bordes_templados(core + "G", region, perdidas3=1) == core
    # una sola base perdida con dos extra: solo se recorta la compensada
    assert _recorte_bordes_templados(
        "GC" + core, region, perdidas5=1) == "C" + core
    # inserción NO templada: no se toca (aunque lleve borde extra)
    azar = "".join(rng.choice("ACGT") for _ in range(20))
    while any(sum(a != b for a, b in zip(azar, region[p:p + 20])) <= 3
              for p in range(181)):
        azar = "".join(rng.choice("ACGT") for _ in range(20))
    assert _recorte_bordes_templados("G" + azar, region) == "G" + azar
    assert _recorte_bordes_templados(
        "G" + azar, region, perdidas5=1) == "G" + azar
    # inserción corta (NHEJ pequeño): no se toca
    assert _recorte_bordes_templados("GATTACAGC", region) == "GATTACAGC"


def test_bordes_de_compensacion_en_insercion_templada(tmp_path):
    """Inserción templada + bases de compensación → un solo grupo.

    Replica el caso real: la inserción de la unión es una copia de la
    región deleción y el flanco 5' lleva un homopolímero (AAA) que el
    ruido de Nanopore lee como AA o A. Esas lecturas pierden 1-2 bases
    en el flanco y el alineador compensa el desplazamiento colocando
    bases fantasma en el borde 5' de la inserción. Tras el recorte de
    bordes templados (basado en el balance de bases perdidas), todas
    deben caer en el mismo grupo mostrando la inserción limpia.
    """
    WT2, DS2, DE2 = 700, 200, 400
    rng = random.Random(21)
    wt = list("".join(rng.choice("ACGT") for _ in range(WT2)))
    wt[190:193] = "AAA"                       # homopolímero del flanco 5'
    wt = "".join(wt)
    _write_fasta(tmp_path / "wt.fa", wt)

    plantilla = wt[330:353]                  # dentro de la deleción
    ins = plantilla[:10] + ("A" if plantilla[10] != "A" else "C") \
        + plantilla[11:]                     # copia con 1 discrepancia
    # base(s) de compensación que NO casa con el final del flanco 5'
    # (para que el alineador las deje dentro de la inserción)
    comp = "G" if wt[DS2 - 2] != "G" else "T"

    # flanco con el homopolímero leído con pérdida: AAA -> AA / A
    flanko_aa = wt[:190] + "AA" + wt[193:DS2 - 1]      # pierde 1 A
    flanko_a = wt[:190] + "A" + wt[193:DS2 - 1]        # pierde 2 A

    reads = ([(f"limpia_{i}", wt[:DS2 - 1] + ins + wt[DE2:])
              for i in range(10)]
             + [(f"comp1_{i}", flanko_aa + comp + ins + wt[DE2:])
                for i in range(6)]
             + [(f"comp2_{i}", flanko_a + comp + "C" + ins + wt[DE2:])
                for i in range(4)]
             + [(f"wt_{i}", wt) for i in range(20)])
    rng.shuffle(reads)
    fq = tmp_path / "muestra4.fastq"
    with open(fq, "w") as f:
        for rid, s in reads:
            f.write(f"@{rid}\n{s}\n+\n{'I' * len(s)}\n")

    res = run_full_pipeline(
        tmp_path / "job4", fq, tmp_path / "wt.fa", "mix4",
        params={
            "deletion_start": DS2, "deletion_end": DE2,
            "min_quality": 10, "min_length": 400, "max_length": 900,
            "breakpoint_tolerance": 20, "min_sv_deletion": 100,
            "nhej_search_window": 200, "wt_min": 650, "wt_max": 760,
        },
        log=lambda m: None,
    )

    delecs = res["deletion_view"]["deleciones"]
    con_ins = [d for d in delecs if d["ins_union"]]
    assert con_ins, "debe haber uniones con inserción"
    # un único grupo con la inserción limpia y las 20 lecturas
    assert len(con_ins) == 1, \
        f"esperaba 1 grupo con inserción, hay {len(con_ins)}: " \
        + "; ".join(f"{d['ins_union'][:30]}(n={d['n']})" for d in con_ins)
    assert con_ins[0]["ins_union"] == ins, con_ins[0]["ins_union"]
    assert con_ins[0]["n"] == 20, con_ins[0]["n"]


def test_bases_reales_en_borde_5_no_se_recortan(tmp_path):
    """Bases extra en el borde 5' sin pérdida de flanco: inserción real.

    Si la lectura no ha perdido bases en un homopolímero del flanco,
    las bases extra en el borde 5' de una inserción templada se
    consideran inserción real: forman su propio grupo y NO se recortan
    ni se fusionan con el grupo de la copia limpia.
    """
    WT2, DS2, DE2 = 700, 200, 400
    rng = random.Random(29)
    wt = list("".join(rng.choice("ACGT") for _ in range(WT2)))
    wt = "".join(wt)
    _write_fasta(tmp_path / "wt.fa", wt)

    plantilla = wt[330:353]
    ins = plantilla[:10] + ("A" if plantilla[10] != "A" else "C") \
        + plantilla[11:]
    # bases reales del borde 5' (sin pérdida en el flanco, no son
    # compensación de homopolímero): con esta semilla el alineador no
    # las absorbe en el flanco ni en la región deleción y quedan dentro
    # del patrón de la inserción
    borde = "GC"

    reads = ([(f"limpia_{i}", wt[:DS2 - 1] + ins + wt[DE2:])
              for i in range(10)]
             + [(f"real_{i}", wt[:DS2 - 1] + borde + ins + wt[DE2:])
                for i in range(6)]
             + [(f"wt_{i}", wt) for i in range(20)])
    rng.shuffle(reads)
    fq = tmp_path / "muestra5.fastq"
    with open(fq, "w") as f:
        for rid, s in reads:
            f.write(f"@{rid}\n{s}\n+\n{'I' * len(s)}\n")

    res = run_full_pipeline(
        tmp_path / "job5", fq, tmp_path / "wt.fa", "mix5",
        params={
            "deletion_start": DS2, "deletion_end": DE2,
            "min_quality": 10, "min_length": 400, "max_length": 900,
            "breakpoint_tolerance": 20, "min_sv_deletion": 100,
            "nhej_search_window": 200, "wt_min": 650, "wt_max": 760,
        },
        log=lambda m: None,
    )

    delecs = res["deletion_view"]["deleciones"]
    con_ins = [d for d in delecs if d["ins_union"]]
    assert len(con_ins) == 2, \
        f"esperaba 2 grupos con inserción, hay {len(con_ins)}: " \
        + "; ".join(f"{d['ins_union'][:30]}(n={d['n']})" for d in con_ins)
    por_ins = {d["ins_union"]: d["n"] for d in con_ins}
    assert por_ins.get(ins) == 10, por_ins
    # las 6 lecturas con las 2 bases 5' reales forman su propio grupo
    assert por_ins.get(borde + ins) == 6, por_ins

