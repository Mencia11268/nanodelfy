"""
Clasificador de lecturas a partir de los dos BAMs (vs WT y vs EDIT).

Clases:
  WT                  → alelo silvestre intacto
  DEL_EXACT           → deleción precisa en los breakpoints esperados (±tol estrecho)
  DEL_EXTENDED        → deleción mayor que la esperada (resección NHEJ extra)
  DEL_TRUNCATED       → deleción menor que la esperada (microhomología)
  DEL_WITH_INSERTION  → deleción + inserción adyacente en la unión (NHEJ con insert)
  INV                 → inversión del fragmento entre breakpoints
  OTHER               → indels imprecisos / NHEJ aberrante / lecturas no informativas

Las cuatro subclases DEL_* se agrupan también como "DEL_ANY" (eficiencia total
de deleción), y todas las ediciones (DEL_* + INV) como "EDITED_ANY".

Parametrización del NHEJ:
  - breakpoint_tolerance (tol_exact): ventana para llamar deleción EXACTA
  - nhej_search_window (search_w):   ventana para reconocer el evento como
                                     deleción imprecisa (típicamente 500–1000 nt
                                     en cada breakpoint). Lo que cae fuera de
                                     esta ventana se considera OTHER.
  - min_junction_insertion: tamaño mínimo de inserción (nt) en la unión para
                            clasificar como DEL_WITH_INSERTION.
"""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple, Dict, Optional
import pandas as pd

from .filter_align import AlnRecord


# Códigos CIGAR de pysam: 0=M 1=I 2=D 3=N 4=S 5=H 6=P 7== 8=X
CONSUMES_REF = {0, 2, 3, 7, 8}
CONSUMES_QUERY_NO_CLIP = {0, 1, 7, 8}


# Etiquetas de clase y orden canónico para los reports
CLASS_ORDER = [
    "WT",
    "DEL_EXACT",
    "DEL_EXTENDED",
    "DEL_TRUNCATED",
    "DEL_WITH_INSERTION",
    "INV",
    "WT_PARTIAL",
    "WT_TRUNCATED",
    "EDIT_TRUNCATED",
    "AMBIGUOUS_TRUNCATED",
    "OTHER",
]
CLASS_DESC = {
    "WT":                  "Alelo silvestre intacto",
    "DEL_EXACT":           "Deleción precisa en los breakpoints esperados",
    "DEL_EXTENDED":        "Deleción extendida (resección adicional por NHEJ)",
    "DEL_TRUNCATED":       "Deleción más corta de lo esperado",
    "DEL_WITH_INSERTION":  "Deleción con inserción en la unión NHEJ",
    "INV":                 "Inversión del fragmento entre breakpoints",
    "WT_PARTIAL":          "Lectura WT que atraviesa al menos un breakpoint sin deleción (informativa)",
    "WT_TRUNCATED":        "Lectura truncada que prefiere WT (no informativa de edición)",
    "EDIT_TRUNCATED":      "Lectura truncada que prefiere EDIT",
    "AMBIGUOUS_TRUNCATED": "Lectura truncada en región idéntica WT/EDIT (no informativa)",
    "OTHER":               "Indels imprecisos / no informativas",
}
DEL_CLASSES = {"DEL_EXACT", "DEL_EXTENDED", "DEL_TRUNCATED", "DEL_WITH_INSERTION"}
EDITED_CLASSES = DEL_CLASSES | {"INV"}
WT_CLASSES = {"WT", "WT_PARTIAL"}
# Las clases TRUNCATED y AMBIGUOUS_TRUNCATED se excluyen del denominador informativo
NON_INFORMATIVE_CLASSES = {"WT_TRUNCATED", "EDIT_TRUNCATED", "AMBIGUOUS_TRUNCATED"}


@dataclass
class Params:
    deletion_start: int            # 1-based, inclusivo
    deletion_end: int              # 1-based, inclusivo
    wt_length: int
    breakpoint_tolerance: int = 20
    min_sv_deletion: int = 1500
    min_edit_ref_coverage: float = 0.90
    min_mapq: int = 20
    # Ventana de longitud de lecturas "tipo WT" (1-based, sobre la lectura):
    # None → se deriva del tamaño del amplicón en __post_init__ (±4 %).
    wt_min: Optional[int] = None
    wt_max: Optional[int] = None
    # NHEJ — nuevos parámetros
    nhej_search_window: int = 500          # ventana de búsqueda alrededor de cada breakpoint
    min_junction_insertion: int = 5        # ≥ 5 nt insertados → DEL_WITH_INSERTION
    min_total_indel_for_junction: int = 3  # I+D adyacentes deben sumar al menos esto
    # Parche R3EB_B: secuencia del amplicón WT (mayúsculas) para reconstruir
    # uniones de deleción a partir de colas soft-clipped. Si es None, esa
    # fuente de evidencia queda desactivada (comportamiento previo).
    wt_seq: Optional[str] = None
    # Capa estadística: calibrar el ruido con las lecturas WT de la muestra
    # (añade p-valores empíricos; no altera la clasificación).
    calibration: bool = True

    def __post_init__(self):
        # Ventana WT: si no se especifica, se deriva del amplicón (±4 %).
        # Para deleciones pequeñas la ventana es orientativa: solo separa
        # lecturas de longitud plenamente WT del resto.
        if self.wt_min is None:
            self.wt_min = round(self.wt_length * 0.96)
        if self.wt_max is None:
            self.wt_max = round(self.wt_length * 1.04)


@dataclass
class ReadEvidence:
    read_id: str
    length: int = 0
    wt_alns: List[AlnRecord] = field(default_factory=list)
    edit_alns: List[AlnRecord] = field(default_factory=list)


# ---------- utilidades CIGAR -----------------------------------------------

def aligned_query_fraction(aln) -> float:
    if not aln.cigartuples:
        return 0.0
    q_aln = sum(l for op, l in aln.cigartuples if op in CONSUMES_QUERY_NO_CLIP)
    q_tot = aln.infer_read_length() or 1
    return q_aln / q_tot


def alignment_score(aln) -> int:
    try:
        return aln.get_tag("AS")
    except KeyError:
        return 0


def best_primary(alns):
    prim = [a for a in alns if not a.is_secondary and not a.is_supplementary]
    return max(prim, key=alignment_score) if prim else None


def find_large_deletions(aln, min_len: int) -> List[Tuple[int, int, int]]:
    """Devuelve [(ref_start, ref_end, length)] para cada D ≥ min_len."""
    out = []
    ref = aln.reference_start
    for op, length in (aln.cigartuples or []):
        if op == 2 and length >= min_len:
            out.append((ref, ref + length, length))
        if op in CONSUMES_REF:
            ref += length
    return out


def find_junction_events(aln, min_d: int) -> List[dict]:
    """
    Busca eventos compuestos en el CIGAR centrados en una D grande.
    Devuelve, para cada D ≥ min_d encontrada:
        {
          'd_start', 'd_end', 'd_len':   coords y longitud de la deleción
          'ins_5p':  inserción inmediatamente 5' de la D (puede ser 0)
          'ins_3p':  inserción inmediatamente 3' de la D
          'ins_total': suma de inserciones adyacentes
        }
    "Adyacente" significa separadas únicamente por operaciones que no consumen
    query/ref (P, H) o por matches cortos de hasta 5 nt — capturamos el caso
    típico del NHEJ donde I y D aparecen pegados o casi pegados en el CIGAR.
    """
    if not aln.cigartuples:
        return []
    cig = aln.cigartuples

    # Pre-calcular posición de referencia y query para cada operación
    ref_positions = []   # ref start de cada op
    qry_positions = []   # query start
    ref = aln.reference_start
    qry = 0
    for op, length in cig:
        ref_positions.append(ref)
        qry_positions.append(qry)
        if op in CONSUMES_REF:
            ref += length
        if op in CONSUMES_QUERY_NO_CLIP or op in (4, 1):
            qry += length

    out = []
    near_match_limit = 5  # un M ≤ 5 nt entre I y D se considera "pegado"

    for i, (op, length) in enumerate(cig):
        if op != 2 or length < min_d:
            continue
        d_start = ref_positions[i]
        d_end = d_start + length

        # Buscar inserción inmediatamente 5' (mirando hacia atrás)
        ins_5p = 0
        j = i - 1
        while j >= 0:
            op_j, len_j = cig[j]
            if op_j == 1:
                ins_5p += len_j
                j -= 1
            elif op_j in (0, 7, 8) and len_j <= near_match_limit:
                j -= 1  # match corto: seguimos mirando
            else:
                break

        # Buscar inserción inmediatamente 3'
        ins_3p = 0
        j = i + 1
        while j < len(cig):
            op_j, len_j = cig[j]
            if op_j == 1:
                ins_3p += len_j
                j += 1
            elif op_j in (0, 7, 8) and len_j <= near_match_limit:
                j += 1
            else:
                break

        out.append({
            "d_start": d_start, "d_end": d_end, "d_len": length,
            "ins_5p": ins_5p, "ins_3p": ins_3p,
            "ins_total": ins_5p + ins_3p,
        })
    return out


# ---------- recogida de evidencias -----------------------------------------

def collect_evidence_from_records(wt_records: List[AlnRecord],
                                  edit_records: List[AlnRecord]) -> Dict[str, ReadEvidence]:
    """Construye el diccionario de evidencias a partir de las listas de AlnRecord
    devueltas por align_minimap2 (no necesita leer BAMs)."""
    ev: Dict[str, ReadEvidence] = {}
    for records, attr in [(wt_records, "wt_alns"), (edit_records, "edit_alns")]:
        for aln in records:
            rec = ev.setdefault(aln.query_name, ReadEvidence(aln.query_name))
            if aln.query_length and not rec.length:
                rec.length = aln.query_length
            if not aln.is_unmapped:
                getattr(rec, attr).append(aln)
    return ev


# ---------- reglas de clasificación ----------------------------------------

def is_inversion(ev: ReadEvidence, p: Params) -> Tuple[bool, dict]:
    """Detecta inversión del fragmento [ds, de] (split-read o aln interno)."""
    ds = p.deletion_start - 1
    de = p.deletion_end
    tol = p.breakpoint_tolerance

    wt = [a for a in ev.wt_alns if not a.is_secondary]

    # ---- Heurística A: split-read canónico
    if len(wt) >= 2:
        forward = [a for a in wt if not a.is_reverse]
        reverse = [a for a in wt if a.is_reverse]
        if forward and reverse:
            fwd_spans = [(a.reference_start, a.reference_end) for a in forward]
            rev_spans = [(a.reference_start, a.reference_end) for a in reverse]
            flanks_ok = any(s < ds + tol for s, _ in fwd_spans) and \
                        any(e > de - tol for _, e in fwd_spans)
            inv_inside = any(
                (s < ds + tol and e > de - tol) or
                (s >= ds - tol and e <= de + tol)
                for s, e in rev_spans
            )
            if flanks_ok and inv_inside:
                return True, {
                    "rule": "INV_split_read",
                    "fwd_segments": len(forward),
                    "rev_segments": len(reverse),
                    "rev_spans": ";".join(f"{s}-{e}" for s, e in rev_spans),
                }

    # ---- Heurística B: read de longitud WT con alineamiento sólo a la región interna
    if p.wt_min <= ev.length <= p.wt_max:
        wt_prim = best_primary(ev.wt_alns)
        edit_prim = best_primary(ev.edit_alns)
        if wt_prim is not None and (wt_prim.mapping_quality or 0) >= p.min_mapq:
            rs, re_ = wt_prim.reference_start, wt_prim.reference_end
            inside_inv_region = (rs >= ds - tol and re_ <= de + tol)
            edit_ref_len = p.wt_length - (de - ds)
            edit_cov = ((edit_prim.reference_length or 0) / max(edit_ref_len, 1)
                        if edit_prim else 0.0)
            if inside_inv_region and edit_cov < p.min_edit_ref_coverage:
                return True, {
                    "rule": "INV_internal_only",
                    "wt_aln_span": f"{rs}-{re_}",
                    "edit_cov": round(edit_cov, 3),
                }

    return False, {}


# ---------- Búsqueda directa de la unión (deleciones pequeñas) -------------
#
# Problema: en deleciones pequeñas (~120 nt) las referencias WT y EDIT son
# casi idénticas, el margen de score EDIT-WT es comparable al ruido de
# secuenciación y el fallback por alineamiento contra EDIT apenas dispara.
# La unión de la deleción (flanco 5' inmediatamente seguido del flanco 3')
# es una señal INDEPENDIENTE del tamaño de la deleción: solo existe en
# lecturas editadas. La buscamos directamente en la secuencia de la lectura
# con semillas exactas + verificación de identidad de los flancos.

_JUNCTION_FLANK = 25      # nt de flanco a cada lado de la unión
_JUNCTION_SEED_K = 8      # longitud de las semillas exactas
_JUNCTION_MH_MAX = 10     # microhomología máxima tolerada (nt)
_JUNCTION_MIN_ID = 0.60   # identidad mínima flanco vs lectura


def _find_junction_in_seq(seq: str, flank5: str, flank3: str, slack: int,
                          k: int = _JUNCTION_SEED_K,
                          flank: int = _JUNCTION_FLANK) -> Optional[int]:
    """Busca la unión (flanco5 seguido de flanco3) en una secuencia.

    Devuelve la separación observada entre flancos (inserción NHEJ, 0 si
    unión exacta, negativo si microhomología) o None si no hay unión.
    """
    # Semillas: dos por flanco, pegadas a la unión.
    # (semilla, distancia desde el inicio de la semilla a la unión)
    #   lado 5': la unión está DESPUÉS del final de la semilla
    #   lado 3': la unión está ANTES del inicio de la semilla
    seeds5 = [(flank5[-2 * k:-k], 2 * k), (flank5[-k:], k)]
    seeds3 = [(flank3[:k], 0), (flank3[k:2 * k], -k)]

    # Índice de k-mers de la secuencia
    idx: Dict[str, List[int]] = {}
    for i in range(0, len(seq) - k + 1):
        idx.setdefault(seq[i:i + k], []).append(i)

    def _identity(a: str, b: str) -> float:
        return sum(x == y for x, y in zip(a, b)) / max(len(b), 1)

    for s5, dist5 in seeds5:
        for hit5 in idx.get(s5, ()):
            junction5 = hit5 + dist5        # posición de la unión (lado 5')
            for s3, dist3 in seeds3:
                for hit3 in idx.get(s3, ()):
                    junction3 = hit3 + dist3  # posición de la unión (lado 3')
                    gap = junction3 - junction5
                    if not (-_JUNCTION_MH_MAX <= gap <= slack):
                        continue
                    # Verificación de identidad de ambos flancos
                    seg5 = seq[junction5 - flank:junction5]
                    seg3 = seq[junction3:junction3 + flank]
                    if len(seg5) < flank or len(seg3) < flank:
                        continue
                    if _identity(seg5, flank5) < _JUNCTION_MIN_ID:
                        continue
                    if _identity(seg3, flank3) < _JUNCTION_MIN_ID:
                        continue
                    return gap
    return None


def _infer_deletion_from_junction_search(ev: ReadEvidence, p: Params,
                                         flank: int = _JUNCTION_FLANK
                                         ) -> Optional[dict]:
    """Cuarta fuente de evidencia: la unión esperada (flanco5+flanco3) está
    literalmente en la secuencia de la lectura.

    Se buscan semillas exactas de 8 nt a ambos lados de la unión y se
    comprueba que la distancia entre flancos sea coherente con la deleción
    esperada (0–slack, tolerando inserciones NHEJ y microhomología corta).
    Un evento aquí tiene breakpoints ASUMIDOS (= esperados), como el
    fallback por alineamiento contra EDIT, pero no depende de que minimap2
    resuelva bien la lectura.

    ESPECIFICIDAD: si la unión ya existe de forma nativa en el amplicón WT
    (p. ej. flancos dentro de repeticiones compuestas como los sitios FRT,
    que son pseudopalíndromos), la señal no es evidencia de edición y la
    búsqueda queda desactivada para todo el amplicón.
    """
    if not p.wt_seq:
        return None
    ds = p.deletion_start - 1
    de = p.deletion_end
    expected_del = de - ds
    k = _JUNCTION_SEED_K
    if expected_del <= 4 * k:      # deleción ínfima: los flancos se solaparían
        return None

    flank5 = p.wt_seq[max(0, ds - flank):ds]
    flank3 = p.wt_seq[de:de + flank]
    if len(flank5) < flank or len(flank3) < flank:
        return None

    # Margen de separación entre flancos: 0 = unión exacta; valores
    # positivos = inserción NHEJ; negativos = microhomología. Debe quedar
    # muy por debajo del tamaño esperado (si no, es una lectura WT normal:
    # sus flancos están a `expected_del` de distancia).
    slack = min(50, expected_del - 4 * k)

    # Comprobación de especificidad (cacheada por Params): ¿existe la
    # unión en el propio WT? En el WT intacto los flancos están a
    # `expected_del`; una aparición con gap pequeño solo puede venir de
    # repeticiones internas → desactivar.
    specific = getattr(p, "_junction_specific", None)
    if specific is None:
        native = _find_junction_in_seq(p.wt_seq, flank5, flank3, slack,
                                      k=k, flank=flank)
        specific = native is None
        p._junction_specific = specific   # type: ignore[attr-defined]
    if not specific:
        return None

    # Secuencia de la lectura (en orientación de referencia) del
    # alineamiento primario; es la misma para WT y EDIT.
    prim = best_primary(ev.wt_alns) or best_primary(ev.edit_alns)
    seq = prim.query_sequence if prim is not None else None
    if not seq:
        return None

    gap = _find_junction_in_seq(seq, flank5, flank3, slack, k=k, flank=flank)
    if gap is None:
        return None
    ins = max(0, gap)
    return {
        "d_start": ds, "d_end": de, "d_len": expected_del,
        "ins_5p": 0, "ins_3p": ins, "ins_total": ins,
        "source": "junction_search",
        "junction_ins": ins,
    }


def _infer_deletion_from_edit_alignment(ev: ReadEvidence, p: Params) -> Optional[dict]:
    """
    Tercera fuente de evidencia: alineamiento contra REF_EDIT.

    Cuando minimap2 no consigue resolver la lectura en dos segmentos contra
    REF_WT (deja un soft-clip grande en un solo alineamiento), sí suele alinear
    la lectura completa contra REF_EDIT (más corta). Si la lectura cubre bien
    la referencia editada y el alineamiento contra WT es muy peor, podemos
    afirmar que se trata de una deleción.

    En ese caso los breakpoints se asumen como los esperados (sin más resolución,
    pero suficientes para clasificar la lectura), y la inserción en la unión se
    estima a partir del exceso de longitud de la query frente al tramo de la
    referencia editada cubierto.
    """
    edit = best_primary(ev.edit_alns)
    if edit is None or (edit.mapping_quality or 0) < p.min_mapq:
        return None

    edit_ref_len = p.wt_length - (p.deletion_end - (p.deletion_start - 1))
    cov = (edit.reference_length or 0) / max(edit_ref_len, 1)
    if cov < p.min_edit_ref_coverage:
        return None

    # Que el alineamiento contra WT sea claramente peor (= AS edit > AS wt)
    wt_prim = best_primary(ev.wt_alns)
    as_edit = alignment_score(edit)
    as_wt = alignment_score(wt_prim) if wt_prim is not None else 0
    if as_edit <= as_wt * 1.1:
        return None  # alineamientos comparables → probablemente no es DEL clara

    ds = p.deletion_start - 1
    de = p.deletion_end
    expected_del = de - ds

    # Estadísticas del alineamiento contra EDIT
    n_insert = sum(l for op, l in (edit.cigartuples or []) if op == 1)
    n_del_total = sum(l for op, l in (edit.cigartuples or []) if op == 2)
    # D "grande" interno = exceso de deleción más allá del esperado
    big_internal_d = sum(l for op, l in (edit.cigartuples or [])
                         if op == 2 and l >= 30)
    # Soft-clips 5' y 3' del alineamiento contra EDIT (indican pérdida de
    # secuencia en los extremos del producto editado)
    cig = edit.cigartuples or []
    soft_5p = cig[0][1] if cig and cig[0][0] == 4 else 0
    soft_3p = cig[-1][1] if cig and cig[-1][0] == 4 else 0

    # Bases "perdidas" de la referencia editada (D grandes + extremos sin cubrir)
    missing_edit = (edit_ref_len - (edit.reference_length or 0)) + big_internal_d
    # Indels netos (insert - del pequeños, sin contar el D grande que ya cuenta
    # como ampliación de la deleción)
    small_indel_net = n_insert - (n_del_total - big_internal_d)

    # Si la deleción fue más extensa (missing_edit > tol): DEL_EXTENDED
    # Si la lectura es más larga que el editado esperado y NO hay missing_edit:
    #    → puede ser DEL_WITH_INSERTION (insert pequeño)
    #    → o DEL_TRUNCATED (insert > ~50 nt, deleción acortada por microhomología)
    # Decidimos d_len y ins por separado:
    # Umbral por defecto: inserciones NHEJ típicas ≤ 100 nt. Por encima de eso
    # se interpreta como deleción acortada (microhomología extendida) más
    # natural que como inserción.
    INS_VS_TRUNC_THRESHOLD = 100

    if missing_edit > p.breakpoint_tolerance:
        d_len_obs = expected_del + missing_edit
        ins_obs = 0
    elif small_indel_net > INS_VS_TRUNC_THRESHOLD:
        # Insert muy grande → más verosímil que la deleción fue más corta
        d_len_obs = expected_del - small_indel_net
        ins_obs = 0
    elif small_indel_net >= p.min_junction_insertion:
        # Insert pequeño-moderado en la unión (1–100 nt) → DEL_WITH_INSERTION
        d_len_obs = expected_del
        ins_obs = small_indel_net
    else:
        d_len_obs = expected_del
        ins_obs = 0

    return {
        "d_start": ds, "d_end": de,
        "d_len": d_len_obs,
        "ins_5p": 0, "ins_3p": ins_obs,
        "ins_total": ins_obs,
        "missing_edit": max(0, missing_edit),
        "big_internal_d": big_internal_d,
        "small_indel_net": small_indel_net,
        "soft_5p": soft_5p, "soft_3p": soft_3p,
        "source": "edit_alignment",
        "edit_cov": round(cov, 3),
    }


# ---------- Parche R3EB_B: uniones reconstruidas desde soft-clips ----------
#
# Problema: cuando la deleción deja un flanco corto en uno de los extremos
# del producto (p. ej. R3EB_B: unión 1145|3068 deja solo ~49 nt tras la
# unión), minimap2 no puede re-alinear ese flanco corto ni como
# alineamiento suplementario (no genera SA). La lectura queda como
# alineamiento primario de la parte 5' con un soft-clip 3', sin ningún D
# grande en el CIGAR, y como la parte alineada es idéntica en WT y EDIT
# (AS_wt == AS_edit) caía en AMBIGUOUS_TRUNCATED.
#
# Solución: si el alineamiento primario contra WT termina (o empieza) cerca
# de un breakpoint esperado y la cola soft-clipped (≥ 25 nt) encaja en la
# referencia más allá del OTRO breakpoint (≥ min_sv_deletion nt de salto),
# reconstruimos la unión: la deleción observada va desde el final del bloque
# alineado hasta la posición donde encaja la cola.

def _kmer_index(ref_seq: str, k: int = 8) -> Dict[str, List[int]]:
    """Índice k-mer de la referencia WT (amplicón pequeño: coste trivial)."""
    from collections import defaultdict
    idx: Dict[str, List[int]] = defaultdict(list)
    for i in range(len(ref_seq) - k + 1):
        idx[ref_seq[i:i + k]].append(i)
    return dict(idx)


def _place_seq_in_ref(seq: str, ref_seq: str, kmer_idx: Dict[str, List[int]],
                      k: int = 8) -> Optional[Tuple[int, int]]:
    """Coloca `seq` (≥ 25 nt) en la referencia mediante índice k-mer +
    distancia de edición (DP). Devuelve (posicion_inicio, errores) del mejor
    emparejamiento de prefijo, o None si no hay emparejamiento fiable.

    Colocar el PREFIJO de la secuencia basta para fijar la unión: en las
    colas 3' soft-clipped la unión está justo al inicio de la cola.
    """
    if len(seq) < 25 or not ref_seq:
        return None
    best = None
    for L in sorted({len(seq), min(len(seq), 45), 30, 25}, reverse=True):
        if L < 25 or L > len(seq):
            continue
        prefix = seq[:L]
        cand: Dict[int, int] = {}
        for off in (2, 8, 15, 22, 30):
            if off + k > L:
                break
            for pos in kmer_idx.get(prefix[off:off + k], []):
                start = pos - off
                if start < 0 or start + L > len(ref_seq):
                    continue
                cand[start] = cand.get(start, 0) + 1
        if not cand:
            continue
        # DP de distancia de edición sobre los 20 mejores candidatos
        for start, votes in sorted(cand.items(), key=lambda x: -x[1])[:20]:
            if votes < 2 and len(cand) > 3:
                continue
            prev = list(range(L + 1))
            for i in range(1, L + 1):
                cur = [i] + [0] * L
                for j in range(1, L + 1):
                    cur[j] = min(prev[j] + 1, cur[j - 1] + 1,
                                 prev[j - 1] + (prefix[i - 1] != ref_seq[start + j - 1]))
                prev = cur
            d = prev[L]
            thr = max(6, int(0.2 * L))
            if d <= thr and (best is None or d < best[1]):
                best = (start, d)
        if best is not None:
            break
    return best


def _infer_deletion_from_softclip(ev: ReadEvidence, p: Params) -> List[dict]:
    """Reconstruye uniones de deleción a partir de soft-clips del alineamiento
    primario contra WT. Devuelve 0, 1 o 2 eventos (cola 3' y/o cabeza 5').

    Requisitos para aceptar la unión:
      - cola/cabeza soft-clipped de ≥ 25 nt;
      - el salto implícito en la referencia ≥ min_sv_deletion nt;
      - (el chequeo de ventana NHEJ alrededor de los breakpoints esperados
         lo aplica classify_deletion, igual que para cigar_D y split_read).
    """
    if not p.wt_seq:
        return []
    wt = best_primary(ev.wt_alns)
    if wt is None or wt.is_unmapped or (wt.mapping_quality or 0) < p.min_mapq:
        return []
    if not wt.cigartuples or wt.query_sequence is None:
        return []
    ref_seq = p.wt_seq
    kmer_idx = getattr(p, "_kmer_idx", None)
    if kmer_idx is None:
        kmer_idx = _kmer_index(ref_seq)
        try:
            object.__setattr__(p, "_kmer_idx", kmer_idx)
        except AttributeError:
            pass

    out: List[dict] = []
    cig = wt.cigartuples

    # --- cola 3' soft-clipped: la lectura sigue en la referencia después
    #     del bloque alineado (salto hacia delante = deleción) ---
    clip3 = cig[-1][1] if cig and cig[-1][0] == 4 else 0
    if clip3 >= 25:
        tail = wt.query_sequence[wt.query_alignment_end:]
        if len(tail) >= 25:
            m = _place_seq_in_ref(tail, ref_seq, kmer_idx)
            if m is not None:
                d_start = wt.reference_end          # fin del bloque alineado
                d_end = m[0]                        # inicio de la cola en ref
                if d_end - d_start >= p.min_sv_deletion:
                    out.append({
                        "d_start": d_start, "d_end": d_end,
                        "d_len": d_end - d_start,
                        "ins_5p": 0, "ins_3p": 0, "ins_total": 0,
                        "source": "softclip_junction",
                        "tail_errors": m[1],
                    })

    # --- cabeza 5' soft-clipped: la lectura empieza ANTES del bloque
    #     alineado y salta hacia delante (deleción entre cabeza y bloque) ---
    clip5 = cig[0][1] if cig and cig[0][0] == 4 else 0
    if clip5 >= 25:
        head = wt.query_sequence[:wt.query_alignment_start]
        if len(head) >= 25:
            # La unión está al FINAL de la cabeza: anclamos exactamente
            # los últimos 25 nt de la cabeza en la referencia.
            suffix = head[-25:]
            m = _place_seq_in_ref(suffix, ref_seq, kmer_idx)
            if m is not None:
                # La cabeza termina en ref m[0]+25; la deleción va desde
                # ahí hasta el inicio del bloque alineado.
                d_start = m[0] + len(suffix)
                d_end = wt.reference_start
                if d_end - d_start >= p.min_sv_deletion:
                    out.append({
                        "d_start": d_start, "d_end": d_end,
                        "d_len": d_end - d_start,
                        "ins_5p": 0, "ins_3p": 0, "ins_total": 0,
                        "source": "softclip_junction",
                        "tail_errors": m[1],
                    })
    return out


def _infer_deletion_from_split(ev: ReadEvidence, p: Params) -> Optional[dict]:
    """
    Cuando minimap2 alinea la lectura en dos segmentos (5' y 3') del mismo
    strand, inferimos la deleción a partir de:
      - gap en la referencia entre el final del 5' y el inicio del 3' = D
      - gap en la query entre los segmentos:
            > 0  → hay una inserción entre los flancos (NHEJ con insert)
            < 0  → microhomología (solapamiento de los flancos)
    Devuelve un dict como find_junction_events() o None.
    """
    alns = [a for a in ev.wt_alns
            if not a.is_secondary and not a.is_unmapped]
    if len(alns) < 2:
        return None

    # Quedarnos con el strand mayoritario para no confundirlo con una inversión
    fwd = [a for a in alns if not a.is_reverse]
    rev = [a for a in alns if a.is_reverse]
    group = fwd if len(fwd) >= len(rev) else rev
    if len(group) < 2:
        return None
    is_rev_group = group[0].is_reverse

    # Para el strand forward, ordenar por query_start; para reverse, al revés.
    # Así left = flanco 5' del read, right = flanco 3' del read.
    group_sorted = sorted(
        group,
        key=lambda a: a.query_alignment_start,
        reverse=is_rev_group,
    )
    left, right = group_sorted[0], group_sorted[-1]
    if left is right:
        return None

    # En la referencia, los segmentos contiguos del read están separados por
    # el hueco de la deleción. Si el strand es reverse, la posición en la
    # referencia del 'left' (5' del read) es la más alta.
    if is_rev_group:
        d_start = right.reference_end
        d_end = left.reference_start
    else:
        d_start = left.reference_end
        d_end = right.reference_start
    if d_end <= d_start:
        return None
    d_len = d_end - d_start

    query_gap = right.query_alignment_start - left.query_alignment_end
    ins_total = max(0, query_gap)

    # --- Ambigüedad por solapamiento en query (repetición en la referencia) ---
    # Si los dos alineamientos reclaman los mismos nt de la lectura
    # (query overlap > 0), la posición exacta de la unión es AMBIGUA: varias
    # posiciones j en [rqs, lqe] producen la MISMA secuencia de lectura
    # (típico en cortes dentro de repeticiones, p. ej. los FRT). Entre las
    # interpretaciones equivalentes elegimos la unión esperada [ds, de)
    # cuando es una de las válidas; si no, se mantiene la observada.
    # (Solo grupo forward: la geometría reverse de una unión estándar no
    # llega a esta inferencia.)
    ds = p.deletion_start - 1
    de = p.deletion_end
    overlap = -query_gap if query_gap < 0 else 0
    snapped = False
    if overlap > 0 and not is_rev_group:
        lqe = left.query_alignment_end
        rqs = right.query_alignment_start

        def _qpos_of_ref(aln, ref_pos):
            """Query position que mapea a ref_pos recorriendo el CIGAR.
            Si ref_pos cae dentro de una D, devuelve el borde izquierdo."""
            ref = aln.reference_start
            q = aln.query_alignment_start
            for op, length in aln.cigartuples:
                if op in (0, 7, 8):  # M, =, X
                    if ref <= ref_pos < ref + length:
                        return q + (ref_pos - ref)
                    ref += length
                    q += length
                elif op == 1:  # I: solo query
                    q += length
                elif op in (2, 3):  # D, N: solo referencia
                    if ref <= ref_pos < ref + length:
                        return q  # borde izquierdo del hueco
                    ref += length
            return None

        def _ref_of_qpos(aln, qpos):
            """Posición de referencia que corresponde a qpos (CIGAR-consciente).
            Si qpos cae dentro de una I, devuelve la posición de referencia
            justo antes de la inserción."""
            ref = aln.reference_start
            q = aln.query_alignment_start
            for op, length in aln.cigartuples:
                if op in (0, 7, 8):  # M
                    if q <= qpos < q + length:
                        return ref + (qpos - q)
                    ref += length
                    q += length
                elif op == 1:  # I
                    if q <= qpos < q + length:
                        return ref  # posición antes de la inserción
                    q += length
                elif op in (2, 3):  # D
                    ref += length
            return None

        # j_exp: posición de query de la unión bajo el alineamiento IZQUIERDO
        # (dónde el alineamiento izquierdo alcanza la referencia ds)
        j_exp = _qpos_of_ref(left, ds)
        if j_exp is not None and rqs <= j_exp < lqe:
            # Verificar que bajo el alineamiento DERECHO, esa misma posición
            # de query corresponde a la referencia de
            ref_at_j = _ref_of_qpos(right, j_exp)
            if ref_at_j == de:
                d_start, d_end = ds, de
                d_len = d_end - d_start
                snapped = True

                # Detectar inserciones EN la unión (dentro del solapamiento):
                # una inserción en la unión aparece como I en el CIGAR de
                # AMBOS alineamientos, en la región de solapamiento. Se
                # cuenta UNA vez (desde el CIGAR izquierdo) para no
                # doble-contar. Esto evita que el snap "borre" la inserción.
                junction_ins = 0
                ref_c = left.reference_start
                q_c = left.query_alignment_start
                for op, length in left.cigartuples:
                    if op in (0, 7, 8):  # M
                        ref_c += length
                        q_c += length
                    elif op == 1:  # I: solo query
                        # inserción dentro del solapamiento [rqs, lqe)
                        if rqs <= q_c and (q_c + length) <= lqe:
                            junction_ins += length
                        q_c += length
                    elif op in (2, 3):  # D
                        ref_c += length
                if junction_ins > 0:
                    ins_total = junction_ins

    return {
        "d_start": d_start, "d_end": d_end, "d_len": d_len,
        "ins_5p": 0, "ins_3p": ins_total,
        "ins_total": ins_total,
        "source": "split_read",
        "query_gap": query_gap,
        "snapped_to_expected": snapped,
    }


def classify_deletion(ev: ReadEvidence, p: Params) -> Tuple[Optional[str], dict]:
    """
    Si la lectura tiene una deleción consistente con la edición esperada,
    devuelve (subclase, info). Si no, devuelve (None, {}).

    Subclases: DEL_EXACT / DEL_EXTENDED / DEL_TRUNCATED / DEL_WITH_INSERTION

    La deleción se infiere de dos fuentes (cualquiera dispara la clasificación):
      1) Evento D grande dentro del CIGAR del alineamiento primario.
      2) Split-read: dos segmentos del mismo strand que flanquean la región
         esperada (caso típico cuando minimap2 parte la lectura en 5' + 3').
    """
    ds = p.deletion_start - 1
    de = p.deletion_end
    expected_len = de - ds
    tol_exact = p.breakpoint_tolerance
    search_w = p.nhej_search_window
    ins_threshold = p.min_junction_insertion

    events: List[dict] = []
    wt = best_primary(ev.wt_alns)
    if wt is not None and (wt.mapping_quality or 0) >= p.min_mapq:
        events.extend(find_junction_events(wt, p.min_sv_deletion))

    split_ev = _infer_deletion_from_split(ev, p)
    if split_ev is not None and split_ev["d_len"] >= p.min_sv_deletion:
        events.append(split_ev)

    # Parche R3EB_B: si no hubo CIGAR-D ni split-read, intentar reconstruir
    # la unión a partir de colas soft-clipped (da breakpoints OBSERVADOS,
    # así que se prefiere al fallback por alineamiento contra EDIT).
    if not events:
        events.extend(_infer_deletion_from_softclip(ev, p))

    # Si todavía no hay eventos, búsqueda directa de la unión en la
    # secuencia de la lectura (crítico en deleciones pequeñas, donde el
    # margen de score EDIT-WT apenas supera el ruido).
    if not events:
        junction_ev = _infer_deletion_from_junction_search(ev, p)
        if junction_ev is not None:
            events.append(junction_ev)

    # Si todavía no hay eventos, recurrimos al alineamiento contra EDIT
    if not events:
        edit_ev = _infer_deletion_from_edit_alignment(ev, p)
        if edit_ev is not None:
            events.append(edit_ev)

    if not events:
        return None, {}

    # Quedarse con el evento cuyo breakpoint está más cerca del esperado
    def dist_to_expected(e):
        return abs(e["d_start"] - ds) + abs(e["d_end"] - de)
    ev_best = min(events, key=dist_to_expected)

    rs, re_, ln = ev_best["d_start"], ev_best["d_end"], ev_best["d_len"]
    bp5_offset = rs - ds
    bp3_offset = re_ - de
    ins_total = ev_best.get("ins_total", 0)
    source = ev_best.get("source", "cigar_D")

    # En cigar_D / split_read / softclip_junction, los breakpoints son
    # observados; chequeamos la ventana NHEJ. Para edit_alignment los
    # breakpoints son asumidos = esperados, así que no aplica.
    if source in ("cigar_D", "split_read", "softclip_junction"):
        if abs(bp5_offset) > search_w or abs(bp3_offset) > search_w:
            return None, {}

    info = {
        "obs_start": rs, "obs_end": re_, "obs_len": ln,
        "bp5_offset": bp5_offset, "bp3_offset": bp3_offset,
        "ins_5p": ev_best.get("ins_5p", 0),
        "ins_3p": ev_best.get("ins_3p", 0),
        "ins_total": ins_total,
        "source": source,
    }
    if source == "edit_alignment":
        info["edit_cov"] = ev_best.get("edit_cov", 0)
        info["missing_edit"] = ev_best.get("missing_edit", 0)
        info["big_internal_d"] = ev_best.get("big_internal_d", 0)

    # ---- Decidir subclase ----
    delta_len = ln - expected_len

    # 1) Inserción en la unión NHEJ (con tamaño de deleción aproximadamente esperado)
    if ins_total >= ins_threshold and abs(delta_len) <= tol_exact:
        info["rule"] = "DEL_with_insertion"
        return "DEL_WITH_INSERTION", info

    # 2) Deleción más extensa de lo esperado (resección NHEJ adicional)
    if delta_len > tol_exact:
        info["rule"] = "DEL_extended"
        return "DEL_EXTENDED", info

    # 3) Deleción más corta de lo esperado (microhomología)
    if delta_len < -tol_exact:
        info["rule"] = "DEL_truncated"
        return "DEL_TRUNCATED", info

    # 4) Deleción exacta
    if source == "cigar_D" and \
       (abs(bp5_offset) > tol_exact or abs(bp3_offset) > tol_exact):
        info["rule"] = "DEL_shifted"  # tamaño ok pero desplazada
        return "DEL_EXACT", info
    info["rule"] = "DEL_exact"
    return "DEL_EXACT", info


def is_wt(ev: ReadEvidence, p: Params) -> bool:
    wt = best_primary(ev.wt_alns)
    if wt is None or (wt.mapping_quality or 0) < p.min_mapq:
        return False
    if find_large_deletions(wt, p.min_sv_deletion):
        return False
    if aligned_query_fraction(wt) < 0.85:
        return False
    if any(a.is_supplementary for a in ev.wt_alns):
        return False
    return p.wt_min <= ev.length <= p.wt_max


def classify_truncated(ev: ReadEvidence, p: Params) -> Tuple[Optional[str], dict]:
    """
    Distingue lecturas TRUNCADAS (que no cubren el amplicón completo) según
    si su mejor alineamiento prefiere WT, EDIT o ambos por igual.

    Esto es crítico en amplicones largos secuenciados con Nanopore, donde un
    porcentaje muy alto de lecturas se trunca por debajo del tamaño esperado.
    Estas lecturas NO permiten distinguir WT de EDIT si no atraviesan el sitio
    de corte, así que las separamos como "no informativas".
    """
    wt = best_primary(ev.wt_alns)
    ed = best_primary(ev.edit_alns)
    as_wt = alignment_score(wt) if wt is not None else 0
    as_ed = alignment_score(ed) if ed is not None else 0

    # Si no hay ningún alineamiento con MAPQ mínima, no es "truncada" sino basura
    if (wt is None or (wt.mapping_quality or 0) < p.min_mapq) and \
       (ed is None or (ed.mapping_quality or 0) < p.min_mapq):
        return None, {}

    # ¿La lectura atraviesa alguno de los breakpoints en el WT?
    ds = p.deletion_start - 1
    de = p.deletion_end
    tol = p.breakpoint_tolerance
    crosses_bp = False
    if wt is not None:
        rs, re_ = wt.reference_start, wt.reference_end
        if (rs < ds - tol and re_ > ds + tol) or (rs < de - tol and re_ > de + tol):
            crosses_bp = True

    # Si atraviesa un breakpoint NO es truncada en el sentido informativo:
    # ya debería haber caído en WT, DEL_* o INV antes de llegar aquí.
    # Si llega aquí cruzando un BP, es un caso raro → OTHER (no truncada).
    if crosses_bp:
        return None, {}

    # Decisión por score
    info = {"AS_wt": as_wt, "AS_edit": as_ed, "length": ev.length}
    margin = 0.05  # 5 % de margen para considerar diferencia real
    if as_ed > as_wt * (1 + margin):
        return "EDIT_TRUNCATED", {**info, "rule": "truncated_prefers_edit"}
    if as_wt > as_ed * (1 + margin):
        return "WT_TRUNCATED", {**info, "rule": "truncated_prefers_wt"}
    return "AMBIGUOUS_TRUNCATED", {**info, "rule": "truncated_identical_region"}


def classify_read(ev: ReadEvidence, p: Params) -> Tuple[str, dict]:
    # 1) Inversión tiene prioridad (estructura compleja)
    ok, info = is_inversion(ev, p)
    if ok:
        return "INV", info
    # 2) Deleción (con subclasificación NHEJ)
    sub, info = classify_deletion(ev, p)
    if sub is not None:
        return sub, info
    # 3) Wild-type completo
    if is_wt(ev, p):
        return "WT", {"rule": "WT_colinear"}
    # 4) WT parcial: atraviesa al menos un breakpoint, prefiere WT y no hay
    #    evidencia de deleción grande. Es informativa (es WT no editado) aunque
    #    no cubra el amplicón entero.
    wt = best_primary(ev.wt_alns)
    ed = best_primary(ev.edit_alns)
    if wt is not None and (wt.mapping_quality or 0) >= p.min_mapq \
       and not find_large_deletions(wt, p.min_sv_deletion):
        ds = p.deletion_start - 1
        de = p.deletion_end
        tol = p.breakpoint_tolerance
        rs, re_ = wt.reference_start, wt.reference_end
        crosses_bp = (rs < ds - tol and re_ > ds + tol) or \
                     (rs < de - tol and re_ > de + tol)
        as_wt = alignment_score(wt)
        as_ed = alignment_score(ed) if ed is not None else 0
        if crosses_bp and as_wt >= as_ed * 1.05:
            return "WT_PARTIAL", {
                "rule": "WT_partial_crosses_bp",
                "AS_wt": as_wt, "AS_edit": as_ed,
                "length": ev.length,
            }
    # 5) Truncada: distinguir si prefiere WT, EDIT o es ambigua
    sub, info = classify_truncated(ev, p)
    if sub is not None:
        return sub, info
    # 6) Resto
    return "OTHER", {}


# ---------- API público ----------------------------------------------------

def run_classification(wt_records, edit_records, params: Params,
                       out_tsv: Path, out_summary: Path,
                       out_xlsx: Optional[Path] = None) -> dict:
    ev = collect_evidence_from_records(wt_records, edit_records)
    rows = []
    counts = Counter()
    for rid, e in ev.items():
        cls, info = classify_read(e, params)
        counts[cls] += 1
        rows.append({
            "read_id": rid,
            "length": e.length,
            "class": cls,
            "n_wt_alns": len(e.wt_alns),
            "n_edit_alns": len(e.edit_alns),
            "AS_wt": alignment_score(best_primary(e.wt_alns)) if best_primary(e.wt_alns) else "",
            "AS_edit": alignment_score(best_primary(e.edit_alns)) if best_primary(e.edit_alns) else "",
            "rule": info.get("rule", ""),
            "obs_start": info.get("obs_start", ""),
            "obs_end": info.get("obs_end", ""),
            "obs_len": info.get("obs_len", ""),
            "bp5_offset": info.get("bp5_offset", ""),
            "bp3_offset": info.get("bp3_offset", ""),
            "ins_total": info.get("ins_total", ""),
            "details": ";".join(
                f"{k}={v}" for k, v in info.items()
                if k not in ("rule", "obs_start", "obs_end", "obs_len",
                             "bp5_offset", "bp3_offset", "ins_total")
            ),
        })

    df = pd.DataFrame(rows)

    # ---- Capa estadística opcional: calibración del ruido con las lecturas
    # WT de la muestra (no altera las clases; añade p-valores empíricos).
    calib_params = getattr(params, "calibration", True)
    calib_summary = None
    if calib_params:
        from .calibration import annotate_with_significance
        classes_by_read = {r["read_id"]: r["class"] for r in rows}
        try:
            calib_summary = annotate_with_significance(rows, ev, classes_by_read)
        except Exception as e:            # la capa nunca rompe el pipeline
            calib_summary = None
            print(f"  [AVISO] calibración estadística desactivada: {e}")
        if calib_summary is None:
            print("  [AVISO] calibración estadística omitida: "
                  f"menos de 20 lecturas WT en la muestra")
        df = pd.DataFrame(rows)

    df.to_csv(out_tsv, sep="\t", index=False)

    total = sum(counts.values()) or 1
    non_info = sum(counts.get(k, 0) for k in NON_INFORMATIVE_CLASSES)
    informative = total - non_info or 1

    rows_summary = [
        {"class": k, "description": CLASS_DESC[k],
         "reads": counts.get(k, 0),
         "frequency_pct": round(100 * counts.get(k, 0) / total, 3),
         "frequency_pct_informative": round(100 * counts.get(k, 0) / informative, 3)
                                       if k not in NON_INFORMATIVE_CLASSES else ""}
        for k in CLASS_ORDER
    ]
    # Agregados
    del_any = sum(counts.get(k, 0) for k in DEL_CLASSES)
    edited_any = sum(counts.get(k, 0) for k in EDITED_CLASSES)
    wt_any = sum(counts.get(k, 0) for k in WT_CLASSES)
    rows_summary.append({
        "class": "WT_ANY",
        "description": "WT completo + WT parcial (sin edición)",
        "reads": wt_any,
        "frequency_pct": round(100 * wt_any / total, 3),
        "frequency_pct_informative": round(100 * wt_any / informative, 3),
    })
    rows_summary.append({
        "class": "DEL_ANY",
        "description": "Cualquier deleción (exacta + imprecisa)",
        "reads": del_any,
        "frequency_pct": round(100 * del_any / total, 3),
        "frequency_pct_informative": round(100 * del_any / informative, 3),
    })
    rows_summary.append({
        "class": "EDITED_ANY",
        "description": "Cualquier edición (deleción + inversión)",
        "reads": edited_any,
        "frequency_pct": round(100 * edited_any / total, 3),
        "frequency_pct_informative": round(100 * edited_any / informative, 3),
    })
    rows_summary.append({
        "class": "INFORMATIVE",
        "description": "Lecturas que cubren el sitio de corte (denominador útil)",
        "reads": informative,
        "frequency_pct": round(100 * informative / total, 3),
        "frequency_pct_informative": 100.0,
    })
    rows_summary.append({
        "class": "TOTAL", "description": "Total de lecturas",
        "reads": total, "frequency_pct": 100.0, "frequency_pct_informative": "",
    })
    if calib_summary is not None:
        rows_summary.append({
            "class": "CALIB_NULL",
            "description": "Lecturas WT usadas para calibrar el ruido (nula empírica)",
            "reads": calib_summary["n_null"],
            "frequency_pct": round(100 * calib_summary["n_null"] / total, 3),
            "frequency_pct_informative": "",
        })
        rows_summary.append({
            "class": "CALIB_SIG",
            "description": f"Lecturas editadas con p < {calib_summary['p_threshold']} "
                           "(significativamente por encima del ruido WT)",
            "reads": calib_summary["n_edited_sig"],
            "frequency_pct": round(100 * calib_summary["n_edited_sig"] / total, 3),
            "frequency_pct_informative": "",
        })
        rows_summary.append({
            "class": "CALIB_EXP_FP",
            "description": "Falsos positivos esperados bajo H0 (suma de p-valores "
                           "de las lecturas editadas)",
            "reads": calib_summary["expected_fp"],
            "frequency_pct": "",
            "frequency_pct_informative": "",
        })
    summary = pd.DataFrame(rows_summary)
    summary.to_csv(out_summary, sep="\t", index=False)

    if out_xlsx is not None:
        with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:
            summary.to_excel(xw, sheet_name="Resumen", index=False)
            df.to_excel(xw, sheet_name="Por_lectura", index=False)

    return {
        "total_reads": total,
        "counts": dict(counts),
        "summary": summary.to_dict(orient="records"),
    }
