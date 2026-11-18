"""Capa estadística de calibración del ruido.

Separa lecturas de señal y de ruido SIN alterar la clasificación
determinista: es una capa aditiva que añade columnas por lectura y una
sección al resumen.

Método (calibración empírica, tipo CRISPRLungo pero no paramétrica):

  1. La hipótesis nula H0 es "la lectura es WT (no editada)". Bajo H0,
     una lectura solo contiene ruido de secuenciación.
  2. La distribución nula se estima con las PROPIAS lecturas WT de la
     muestra (clases WT y WT_PARTIAL): son lecturas reales con el mismo
     ruido de secuenciación que el resto de la muestra. Si la muestra no
     tiene lecturas WT suficientes (mínimo 20), la capa se desactiva con
     un aviso (o se indica una muestra control WT con `wt_control`).
  3. Estadísticos por lectura (continuos, no censurados):
       - max_d_wt:  la mayor operación D del CIGAR del alineamiento
                    primario contra WT (evidencia estructural de deleción)
       - as_margin: AS(EDIT) − AS(WT) de los alineamientos primarios
                    (evidencia por preferencia de alineamiento)
  4. Cada lectura recibe p-valores empíricos un lado
        P(estadístico de una lectura WT ≥ el de esta lectura)
     y un indicador de significación (p < 0.05).

El resumen incluye el "número esperado de falsos positivos": la suma de
los p-valores de las lecturas clasificadas como editadas, es decir, cuántas
de ellas serían compatibles con puro ruido WT.

La clasificación determinista (reglas del clasificador) NO cambia: esta
capa solo la anota con significación estadística.
"""
from __future__ import annotations

from bisect import bisect_left
from typing import Dict, List, Optional

from .classify import (
    Params, ReadEvidence, best_primary, alignment_score,
    WT_CLASSES, EDITED_CLASSES, NON_INFORMATIVE_CLASSES,
)

# Mínimo de lecturas WT para estimar la nula; por debajo, se desactiva.
MIN_NULL_READS = 20
# Umbral de significación por lectura
P_SIG = 0.05


def read_noise_metrics(ev: ReadEvidence) -> Dict[str, int]:
    """Estadísticos de ruido por lectura, independientes de su clase.

    - max_d_wt:  mayor D del CIGAR del alineamiento primario WT (nt, 0 si no hay)
    - as_margin: AS(EDIT) − AS(WT) de los mejores alineamientos primarios
    """
    wt = best_primary(ev.wt_alns)
    ed = best_primary(ev.edit_alns)
    max_d = 0
    if wt is not None and wt.cigartuples:
        max_d = max((L for op, L in wt.cigartuples if op == 2), default=0)
    as_wt = alignment_score(wt) if wt is not None else 0
    as_ed = alignment_score(ed) if ed is not None else 0
    return {"max_d_wt": max_d, "as_margin": as_ed - as_wt}


def build_null_distribution(evidence: Dict[str, ReadEvidence],
                            classes: Dict[str, str]) -> Dict:
    """Distribución nula empírica a partir de las lecturas clasificadas WT.

    Devuelve {'n': int, 'max_d_wt': [ordenados], 'as_margin': [ordenados]}.
    """
    null_d: List[int] = []
    null_m: List[int] = []
    for rid, ev in evidence.items():
        if classes.get(rid) not in WT_CLASSES:
            continue
        m = read_noise_metrics(ev)
        null_d.append(m["max_d_wt"])
        null_m.append(m["as_margin"])
    null_d.sort()
    null_m.sort()
    return {"n": len(null_d), "max_d_wt": null_d, "as_margin": null_m}


def _empirical_p(sorted_vals: List[int], x: int) -> Optional[float]:
    """p-valor empírico un lado: P(nula ≥ x). Con corrección (1+n)."""
    if not sorted_vals:
        return None
    n_ge = len(sorted_vals) - bisect_left(sorted_vals, x)
    return (1 + n_ge) / (1 + len(sorted_vals))


def annotate_with_significance(rows: List[dict],
                               evidence: Dict[str, ReadEvidence],
                               classes: Dict[str, str]) -> Optional[dict]:
    """Añade a cada fila las columnas de calibración estadística.

    Devuelve un dict de resumen de la calibración (o None si no hay nula).
    `rows` se modifica in place: cada dict gana las claves
    max_d_wt, as_margin, p_as_margin, p_max_d, sig_noise.
    """
    null = build_null_distribution(evidence, classes)
    if null["n"] < MIN_NULL_READS:
        return None

    n_edited_sig = 0
    expected_fp = 0.0
    for row in rows:
        ev = evidence.get(row["read_id"])
        if ev is None:
            continue
        # Las lecturas TRUNCADAS (regiones idénticas WT/EDIT) tienen margen ~0
        # por construcción: el p-valor no es comparable y se dejan vacías.
        if classes.get(row["read_id"]) in NON_INFORMATIVE_CLASSES:
            row["max_d_wt"] = ""
            row["as_margin"] = ""
            row["p_as_margin"] = ""
            row["p_max_d"] = ""
            row["sig_noise"] = False
            continue
        m = read_noise_metrics(ev)
        p_m = _empirical_p(null["as_margin"], m["as_margin"])
        p_d = _empirical_p(null["max_d_wt"], m["max_d_wt"])
        row["max_d_wt"] = m["max_d_wt"]
        row["as_margin"] = m["as_margin"]
        row["p_as_margin"] = round(p_m, 5) if p_m is not None else ""
        row["p_max_d"] = round(p_d, 5) if p_d is not None else ""
        row["sig_noise"] = bool(p_m is not None and p_m < P_SIG)
        if classes.get(row["read_id"]) in EDITED_CLASSES:
            if row["sig_noise"]:
                n_edited_sig += 1
            if p_m is not None:
                expected_fp += p_m

    return {
        "n_null": null["n"],
        "n_edited_sig": n_edited_sig,
        "expected_fp": round(expected_fp, 2),
        "p_threshold": P_SIG,
    }
