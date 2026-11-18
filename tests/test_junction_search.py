"""Tests de la búsqueda directa de unión (deleciones pequeñas) y de la
capa de calibración estadística del ruido.

La búsqueda de unión es la solución al fallback EDIT débil en deleciones
pequeñas: la unión flanco5+flanco3 solo existe en lecturas editadas,
sea cual sea el tamaño de la deleción.
"""
import random
from dataclasses import field

from nanodelfy.core.classify import (
    Params, ReadEvidence, classify_read, _infer_deletion_from_junction_search,
)
from nanodelfy.core.filter_align import AlnRecord
from nanodelfy.core.calibration import (
    read_noise_metrics, build_null_distribution, annotate_with_significance,
)


def _mk_record(rid: str, seq: str) -> AlnRecord:
    """AlnRecord sintético: alineamiento completo tipo match."""
    return AlnRecord(
        query_name=rid, flag=0, reference_name="amp",
        reference_start=0, mapping_quality=60,
        cigartuples=[(0, len(seq))], query_length=len(seq),
        query_sequence=seq,
    )


def _mk_params(wt: str, del_start: int, del_end: int) -> Params:
    return Params(deletion_start=del_start, deletion_end=del_end,
                  wt_length=len(wt), wt_seq=wt,
                  min_sv_deletion=100, breakpoint_tolerance=20,
                  nhej_search_window=200)


def _mutate(seq: str, rate: float, rng: random.Random) -> str:
    out = []
    for ch in seq:
        r = rng.random()
        if r < rate:                     # sustitución
            out.append(rng.choice("ACGT".replace(ch, "")))
        else:
            out.append(ch)
    return "".join(out)


def _mk_ev(rid: str, seq: str) -> ReadEvidence:
    ev = ReadEvidence(rid, length=len(seq))
    ev.wt_alns = [_mk_record(rid, seq)]
    ev.edit_alns = [_mk_record(rid, seq)]
    return ev


# ---------------------------------------------------------------- unión ---
def test_junction_detectada_en_deleccion_pequena():
    """Lectura EDIT con la deleción exacta de 120 nt: la unión se encuentra
    aunque el CIGAR no tenga ningún D grande (fallback de minimap2)."""
    rng = random.Random(7)
    wt = "".join(rng.choice("ACGT") for _ in range(2000))
    ds, de = 1000, 1120                       # deleción de 120 nt
    edit_read = wt[:1000] + wt[1120:]         # sin la región deletada

    p = _mk_params(wt, ds + 1, de)
    ev = _mk_ev("edit_1", edit_read)
    ev.wt_alns[0].cigartuples = [(0, 1000), (1, 120), (0, 880)]  # I en WT
    # AS del alineamiento contra EDIT > AS contra WT (simulado con tags)
    ev.wt_alns[0]._tags = {}
    ev.edit_alns[0]._tags = {}

    hit = _infer_deletion_from_junction_search(ev, p)
    assert hit is not None
    assert hit["source"] == "junction_search"
    assert hit["d_len"] == 120
    assert hit["ins_total"] == 0


def test_junction_con_ruido_de_secuenciacion():
    """La unión se encuentra con ~3% de errores en la lectura."""
    rng = random.Random(11)
    wt = "".join(rng.choice("ACGT") for _ in range(2000))
    ds, de = 1000, 1120
    edit_read = _mutate(wt[:1000] + wt[1120:], 0.03, rng)
    p = _mk_params(wt, ds + 1, de)
    ev = _mk_ev("edit_ruidosa", edit_read)
    hit = _infer_deletion_from_junction_search(ev, p)
    assert hit is not None


def test_junction_con_insercion_nhej():
    """Unión con 30 nt insertados (NHEJ): detectada, ins=30."""
    rng = random.Random(13)
    wt = "".join(rng.choice("ACGT") for _ in range(2000))
    ds, de = 1000, 1120
    ins = "".join(rng.choice("ACGT") for _ in range(30))
    edit_read = wt[:1000] + ins + wt[1120:]
    p = _mk_params(wt, ds + 1, de)
    ev = _mk_ev("edit_nhej", edit_read)
    hit = _infer_deletion_from_junction_search(ev, p)
    assert hit is not None
    assert hit["ins_total"] == 30


def test_lectura_wt_no_da_union():
    """Lectura WT completa: los flancos están a 120 nt → no hay unión."""
    rng = random.Random(17)
    wt = "".join(rng.choice("ACGT") for _ in range(2000))
    ds, de = 1000, 1120
    p = _mk_params(wt, ds + 1, de)
    ev = _mk_ev("wt_1", wt)
    assert _infer_deletion_from_junction_search(ev, p) is None


def test_wt_con_ruido_no_da_union():
    """Lectura WT con errores: sigue sin detectarse la unión (0 FP)."""
    rng = random.Random(19)
    wt = "".join(rng.choice("ACGT") for _ in range(2000))
    ds, de = 1000, 1120
    noisy_wt = _mutate(wt, 0.03, rng)
    p = _mk_params(wt, ds + 1, de)
    ev = _mk_ev("wt_ruidosa", noisy_wt)
    assert _infer_deletion_from_junction_search(ev, p) is None


def test_union_no_especifica_en_repeticiones_frt():
    """Amplicón con sitios FRT (repetición compuesta): la unión existe de
    forma nativa en el WT → la búsqueda queda desactivada (0 FP, caso
    real MmC7 deleción 1438-3295 con 8.579 falsos positivos sin este
    control de especificidad).
    """
    # Fragmento FRT (pseudopalíndromo compuesto, como el real de MmC7)
    frt = "GAAGTTCCTATTCCGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCGG"
    # La unión "CTAGAAAG|TATAGGAA" existe dentro del propio FRT;
    # los breakpoints se colocan justo ahí, como en el MmC7 real.
    cut = frt.find("CTAGAAAG") + len("CTAGAAAG")     # posición de la unión
    rng = random.Random(31)
    rnd = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    frt1_start, spacer = 500, 1200
    wt = rnd(frt1_start) + frt + rnd(spacer) + frt + rnd(500)
    ds = frt1_start + cut
    de = frt1_start + len(frt) + spacer + cut
    assert ds < de and wt[ds-8:ds] == "CTAGAAAG" and wt[de:de+8] == "TATAGGAA"

    p = _mk_params(wt, ds + 1, de)
    ev = _mk_ev("wt_frt", wt)
    assert _infer_deletion_from_junction_search(ev, p) is None
    # Y también para una lectura WT con ruido
    ev2 = _mk_ev("wt_frt_ruidosa", _mutate(wt, 0.02, rng))
    assert _infer_deletion_from_junction_search(ev2, p) is None


# --------------------------------------------------------- calibración ---
def test_calibracion_separa_senal_de_ruido():
    """Las lecturas WT dan p-valores altos y las EDIT bajos (nula empírica)."""
    rng = random.Random(23)
    wt = "".join(rng.choice("ACGT") for _ in range(2000))
    ds, de = 1000, 1120

    evidence, classes = {}, {}
    # 30 lecturas WT (con ruido) → nula
    for i in range(30):
        seq = _mutate(wt, 0.02, rng)
        ev = _mk_ev(f"wt_{i}", seq)
        ev.wt_alns[0].cigartuples = [(0, len(seq))]
        evidence[f"wt_{i}"] = ev
        classes[f"wt_{i}"] = "WT"
    # 30 lecturas EDIT → señal
    for i in range(30):
        seq = _mutate(wt[:1000] + wt[1120:], 0.02, rng)
        ev = _mk_ev(f"edit_{i}", seq)
        # contra WT: la deleción aparece como I de 120 (peor AS);
        # contra EDIT: alineamiento perfecto (mejor AS)
        ev.wt_alns[0].cigartuples = [(0, 1000), (1, 120), (0, 880)]
        ev.wt_alns[0]._tags["AS"] = 3000
        ev.edit_alns[0]._tags["AS"] = 3500
        evidence[f"edit_{i}"] = ev
        classes[f"edit_{i}"] = "DEL_EXACT"

    rows = [{"read_id": rid, "class": cls} for rid, cls in classes.items()]
    cal = annotate_with_significance(rows, evidence, classes)
    assert cal is not None
    assert cal["n_null"] == 30

    by_id = {r["read_id"]: r for r in rows}
    # WT: p-valor alto (el margen de la WT es ~0, dentro de la nula)
    assert by_id["wt_0"]["p_as_margin"] > 0.5
    # EDIT: margen 500 muy por encima de la nula → p mínimo
    assert by_id["edit_0"]["p_as_margin"] < 0.05
    assert by_id["edit_0"]["sig_noise"] is True
    # FP esperados: suma de p de las 30 EDIT ≈ 30/(30+1)
    assert cal["expected_fp"] < 1.5


def test_calibracion_desactivada_sin_wt():
    """Sin lecturas WT suficientes, la capa se omite (no rompe nada)."""
    rng = random.Random(29)
    wt = "".join(rng.choice("ACGT") for _ in range(500))
    evidence, classes = {}, {}
    for i in range(5):                      # solo 5 WT (< 20)
        ev = _mk_ev(f"wt_{i}", wt)
        evidence[f"wt_{i}"] = ev
        classes[f"wt_{i}"] = "WT"
    rows = [{"read_id": rid, "class": cls} for rid, cls in classes.items()]
    assert annotate_with_significance(rows, evidence, classes) is None
