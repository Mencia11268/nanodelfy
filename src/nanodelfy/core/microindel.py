"""
Detección de microindels (1–50 nt) alrededor de los cortes de guías Cas9.

Añadido a NanoDel como resultado INDEPENDIENTE de la clasificación de grandes
deleciones: una lectura etiquetada WT estructural puede portar un microindel
en el sitio de corte.

Funcionamiento:
  1. Localiza cada guía (20 nt, ARN o ADN) en el amplicón WT, en ambas hebras
     y en todas sus copias (guías contra secuencias repetidas, p.ej. FRT).
     Cas9 corta 3 nt aguas arriba del PAM (NGG): corte romo entre las
     posiciones 17|18 del protospacer.
  2. Alinea las lecturas (mappy, preset map-ont) contra el amplicón WT.
     Para cada sitio de corte se usa el mejor alineamiento (primario o
     suplementario) que ATRAVIESE la ventana del sitio.
  3. Para cada sitio y lectura clasifica:
         NO_INDEL     → ventana limpia (sitio evaluable, sin microindel)
         MICROINDEL   → ≥1 evento I/D de 1 a `max_indel` nt que SOLAPA la
                        ventana (solapamiento por intervalos, no solo por
                        posición inicial)
         OTRO_INDEL   → evento de tamaño intermedio (max_indel < tamaño <
                        large_del_threshold) que solapa la ventana: edición
                        real pero fuera de la definición de microindel
         LARGE_DEL    → el sitio está afectado por una deleción estructural
                        (≥ `large_del_threshold` nt solapando la ventana;
                        incluye la deleción esperada FRT–FRT)
         NOT_EVALUABLE→ la lectura no atraviesa la ventana completa
  4. Agrupa lecturas en ALELOS por los eventos normalizados (representación
     canónica desplazada a la izquierda en contextos repetidos) y extrae el
     contexto de unión RECONSTRUIDO desde la referencia (±30 nt) — etiquetado
     como reconstrucción, no como observación directa de la lectura.

Advertencias de interpretación:
  - Los indels de 1 nt incluyen ruido de secuenciación Nanopore; se reportan
    sin filtrar junto con la calidad media de los eventos.
  - Los porcentajes son condicionados a ventanas localmente evaluables: no
    representan la edición total del locus ni se suman entre sitios.
  - El soporte recurrente (≥3 lecturas) es un análisis descriptivo de patrones;
    no corrige el ruido ni prueba independencia molecular (lecturas de PCR).
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple, Dict, Optional

import mappy as mp

# Códigos CIGAR: 0=M 1=I 2=D 3=N 4=S 5=H 6=P 7== 8=X


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTUacgtu", "TGCATgcata"))[::-1]


def _norm_guide(g: str) -> str:
    """Acepta ARN (U) o ADN (T); devuelve el protospacer en ADN mayúsculas."""
    return g.strip().upper().replace("U", "T")


# ---------------------------------------------------------------- sitios de corte

@dataclass
class CutSite:
    guide_label: str
    guide_seq: str            # protospacer (hebra +)
    strand: str               # '+' o '-'
    protospacer_start: int    # 0-based en la referencia (hebra +)
    pam: Optional[str]        # PAM en la hebra de la guía (NGG cas9 / TTTV cas12a)
    pam_found: bool
    cut_pos: int              # 0-based: corte ENTRE las bases cut_pos-1 y cut_pos.
                              # Para cortes escalonados (cas12a) es el CENTRO de la
                              # región de corte (véase cut_region).
    nuclease: str = "cas9"
    cut_nick_a: Optional[int] = None  # nick PAM-proximal (cortes escalonados)
    cut_nick_b: Optional[int] = None  # nick PAM-distal

    @property
    def cut_region(self) -> Tuple[int, int]:
        """Región de corte (nicks extremos, 0-based 'entre bases').

        cas9  : corte romo → un único punto (ambos valores iguales).
        cas12a: corte escalonado con extremos 5' salientes de 5 nt →
                nicks a 18 y 23 nt del extremo PAM-proximal del protospacer.
        """
        if self.cut_nick_a is None or self.cut_nick_b is None:
            return (self.cut_pos, self.cut_pos)
        lo, hi = sorted((self.cut_nick_a, self.cut_nick_b))
        return (lo, hi)

    @property
    def cut_1based(self) -> str:
        lo, hi = self.cut_region
        if hi == lo:
            return f"entre {lo} y {lo+1}"
        return f"entre {lo+1} y {hi} (región de {hi-lo} nt)"

    @property
    def site_id(self) -> str:
        return f"{self.guide_label}@{self.cut_pos}"


# Parámetros por nucleasa: (longitud esperada del protospacer, tamaño del PAM,
# patrón del PAM en la hebra de la guía, nick PAM-proximal, nick PAM-distal).
# Los nicks se miden en nt desde el extremo PAM-PROXIMAL del protospacer
# (convención "entre bases": un nick a distancia d cae entre las bases
# d-1|d del protospacer, 0-based p+d-1|p+d en la referencia).
_NUCLEASE_SPECS: Dict[str, dict] = {
    # Cas9: corte romo 3 nt aguas arriba del PAM (entre 17|18 del protospacer
    # de 20 nt); PAM NGG inmediatamente 3' del protospacer.
    "cas9": {
        "len_expected": (20,), "pam_len": 3, "pam_pattern": r"^[ACGT]GG$",
        "nicks": (17, 17),
    },
    # Cas12a (Cpf1): corte escalonado con extremos 5' salientes de ~5 nt,
    # nicks a 18 y 23 nt del PAM (Zetsche et al. 2015); PAM TTTV (4 nt)
    # inmediatamente 5' del protospacer; protospacer de 23-24 nt.
    "cas12a": {
        "len_expected": (23, 24), "pam_len": 4, "pam_pattern": r"^TTT[ACGT]$",
        "nicks": (18, 23),
    },
}


def _normalize_nuclease(nuclease: str) -> str:
    n = (nuclease or "cas9").strip().lower().replace("-", "").replace(" ", "")
    aliases = {"cas9": "cas9", "spcas9": "cas9", "cas12a": "cas12a",
               "cpf1": "cas12a", "ascas12a": "cas12a", "lbcas12a": "cas12a"}
    if n not in aliases:
        raise ValueError(f"Nucleasa no soportada: {nuclease!r} "
                         f"(opciones: {', '.join(sorted(set(aliases.values())))})")
    return aliases[n]


def locate_cut_sites(ref_seq: str, guides: Dict[str, str],
                    nuclease: str = "cas9",
                    log=print) -> List[CutSite]:
    """Localiza cada guía en el amplicón (ambas hebras, todas las copias).

    `nuclease` ("cas9" por defecto, o "cas12a") controla la regla de corte y
    el PAM esperado:
      - cas9  : PAM NGG aguas abajo (3'), corte romo 3 nt aguas arriba del PAM.
      - cas12a: PAM TTTV aguas arriba (5'), corte escalonado con nicks a 18 y
                23 nt del extremo PAM-proximal (extremos 5' salientes de 5 nt).
    Para cortes escalonados, `cut_pos` es el centro de la región de corte
    (la ventana de análisis ±window la cubre por completo) y `cut_region`
    devuelve los dos nicks extremos.
    """
    nuc = _normalize_nuclease(nuclease)
    spec = _NUCLEASE_SPECS[nuc]
    L_g, pam_len, pam_pat = spec["len_expected"], spec["pam_len"], spec["pam_pattern"]
    nick_prox, nick_dist = spec["nicks"]
    import re as _re
    pam_re = _re.compile(pam_pat)

    ref = ref_seq.upper()
    sites: List[CutSite] = []
    for label, g_raw in guides.items():
        g = _norm_guide(g_raw)
        if len(g) not in L_g:
            log(f"  [AVISO] {label}: longitud {len(g)} nt "
                f"(se esperaban {' o '.join(map(str, L_g))} para {nuc}); "
                f"se aplican las distancias de corte de {nuc} igualmente")
        # ------- hebra '+': PAM aguas abajo (cas9) o aguas arriba (cas12a)
        start = 0
        while True:
            p = ref.find(g, start)
            if p < 0:
                break
            if nuc == "cas9":
                pam = ref[p + len(g): p + len(g) + pam_len]
            else:  # cas12a: PAM inmediatamente 5' del protospacer
                pam = ref[max(0, p - pam_len): p]
            pam_ok = len(pam) == pam_len and bool(pam_re.match(pam))
            if nuc == "cas9":
                cut = p + nick_prox
                sites.append(CutSite(label, g, "+", p, pam, pam_ok, cut,
                                     nuclease=nuc))
            else:
                nick1, nick2 = p + nick_prox, p + nick_dist
                sites.append(CutSite(label, g, "+", p, pam, pam_ok,
                                    (nick1 + nick2) // 2, nuclease=nuc,
                                    cut_nick_a=nick1, cut_nick_b=nick2))
            start = p + 1
        # ------- hebra '-': espejo (PAM en el extremo opuesto)
        gr = revcomp(g)
        start = 0
        while True:
            p = ref.find(gr, start)
            if p < 0:
                break
            if nuc == "cas9":
                pam = revcomp(ref[max(0, p - pam_len):p]) if p >= pam_len else ""
            else:  # cas12a: PAM aguas abajo en la hebra + (revcomp = TTTV)
                pam = revcomp(ref[p + len(g): p + len(g) + pam_len])
            pam_ok = len(pam) == pam_len and bool(pam_re.match(pam))
            if nuc == "cas9":
                cut = p + len(g) - nick_prox
                sites.append(CutSite(label, g, "-", p, pam, pam_ok, cut,
                                     nuclease=nuc))
            else:
                # nicks a nick_prox/nick_dist del extremo PAM-proximal
                # (el extremo 3' del protospacer en coords de la hebra +)
                nick1 = p + len(g) - nick_prox
                nick2 = p + len(g) - nick_dist
                sites.append(CutSite(label, g, "-", p, pam, pam_ok,
                                    (nick1 + nick2) // 2, nuclease=nuc,
                                    cut_nick_a=nick1, cut_nick_b=nick2))
            start = p + 1
    if not sites:
        raise ValueError("Ninguna de las guías se localiza en el amplicón WT "
                         "(búsqueda exacta en ambas hebras).")
    for s in sites:
        estado = f"PAM {s.pam}" if s.pam_found else \
            f"SIN PAM adyacente (leído: {s.pam or '—'})"
        extra = "" if s.nuclease == "cas9" else f" [{s.nuclease}, escalonado]"
        log(f"  · {s.site_id}: hebra {s.strand}, protospacer "
            f"{s.protospacer_start + 1}-{s.protospacer_start + len(s.guide_seq)}"
            f" (1-based), {estado}, corte {s.cut_1based}{extra}")
    return sites


# ---------------------------------------------------------------- eventos

@dataclass
class IndelEvent:
    op: str                   # 'D' o 'I'
    ref_start: int            # 0-based: D → primera base ref borrada;
                              #         I → inserción justo antes de esta base
    length: int
    ins_seq: str = ""
    mean_q: Optional[float] = None

    def __str__(self):
        if self.op == "D":
            return f"D{self.length}@{self.ref_start + 1}-{self.ref_start + self.length}"
        return f"I{self.length}@antes_de_{self.ref_start + 1}({self.ins_seq})"


def _oriented_fragment(aln, seq: str, qual: Optional[str]):
    """
    Fragmento de lectura en la orientación del alineamiento:
      - hebra + : seq[q_st:q_en]
      - hebra − : revcomp de ese fragmento (y calidad invertida)
    El CIGAR de mappy empieza en q_st (sin soft clips), por lo que el
    recorrido del CIGAR debe hacerse sobre este fragmento desde qpos=0.
    """
    if aln.strand == 1:
        frag = seq[aln.q_st:aln.q_en]
        fqual = qual[aln.q_st:aln.q_en] if qual else None
    else:
        frag = revcomp(seq[aln.q_st:aln.q_en])
        fqual = qual[aln.q_st:aln.q_en][::-1] if qual else None
    return frag, fqual


def _extract_events(aln, frag: str, fqual: Optional[str]):
    """Eventos I/D con coordenadas ref 0-based, sobre el fragmento orientado."""
    events: List[IndelEvent] = []
    m_blocks: List[Tuple[int, int]] = []
    ref = aln.r_st
    qpos = 0
    for length, op in aln.cigar:
        if op in (0, 7, 8):
            m_blocks.append((ref, ref + length))
            ref += length
            qpos += length
        elif op == 1:
            mq = None
            if fqual:
                qs = [ord(c) - 33 for c in fqual[qpos:qpos + length]]
                mq = sum(qs) / len(qs) if qs else None
            events.append(IndelEvent("I", ref, length,
                                     frag[qpos:qpos + length], mq))
            qpos += length
        elif op in (2, 3):
            events.append(IndelEvent("D", ref, length))
            ref += length
    return events, m_blocks


def _deletion_flank_q(fqual: Optional[str], ev_start_in_frag: int,
                      frag_len: int, n: int = 2) -> Optional[float]:
    """Calidad media de ~n bases del fragmento alineado junto a la deleción."""
    if not fqual:
        return None
    lo = max(0, ev_start_in_frag - n // 2)
    hi = min(frag_len, ev_start_in_frag + (n + 1) // 2)
    qs = [ord(c) - 33 for c in fqual[lo:hi]]
    return sum(qs) / len(qs) if qs else None


def _locate_in_fragment(aln, ref_pos: int) -> Optional[int]:
    """Posición (0-based) de ref_pos dentro del fragmento alineado."""
    ref = aln.r_st
    qpos = 0
    for length, op in aln.cigar:
        if op in (0, 7, 8):
            if ref <= ref_pos < ref + length:
                return qpos + (ref_pos - ref)
            ref += length
            qpos += length
        elif op == 1:
            # la inserción consume query (fragmento) pero no referencia
            qpos += length
        elif op in (2, 3):
            if ref <= ref_pos < ref + length:
                return qpos   # cae dentro de la D: borde izquierdo
            ref += length
    return None


# ---------------------------------------------------------------- normalización

def normalize_events(events: List[IndelEvent], ref_seq: str) -> List[IndelEvent]:
    """
    Representación canónica (lo más a la izquierda posible) de cada evento en
    contextos repetidos, para que indels equivalentes agrupen en el mismo alelo.
    """
    out = []
    for ev in events:
        if ev.op == "D":
            s, L = ev.ref_start, ev.length
            while s > 0 and ref_seq[s - 1] == ref_seq[s + L - 1]:
                s -= 1
            out.append(IndelEvent("D", s, L, "", ev.mean_q))
        else:  # I
            pos, S = ev.ref_start, ev.ins_seq
            while pos > 0 and S and S[-1] == ref_seq[pos - 1]:
                S = S[-1] + S[:-1]
                pos -= 1
            out.append(IndelEvent("I", pos, ev.length, S, ev.mean_q))
    return out


# ---------------------------------------------------------------- clasificación

CATEGORIES = ["NO_INDEL", "MICROINDEL", "OTRO_INDEL", "LARGE_DEL",
              "NOT_EVALUABLE"]


@dataclass
class SiteCall:
    site: CutSite
    category: str
    events: List[IndelEvent] = field(default_factory=list)
    large_events: List[IndelEvent] = field(default_factory=list)
    reason: str = ""


def call_site(site: CutSite, aln, events: List[IndelEvent],
              m_blocks: List[Tuple[int, int]], fqual: Optional[str],
              frag_len: int, window: int, large_del_threshold: int,
              min_mapq: int, max_indel: int) -> SiteCall:
    """
    Ventana = intervalo CERRADO [lo, hi] con lo = cut-window, hi = cut+window
    (ambos extremos incluidos):
      - una deleción se considera dentro si SOLAPA el intervalo:
        empieza en una posición <= hi y termina (start+length) > lo
        (puede empezar antes de la ventana y entrar en ella, y puede
        empezar exactamente en hi);
      - una inserción en la posición p (entre las bases p-1 y p) está
        dentro si lo <= p <= hi.
    Si en la ventana hay un microindel (<= max_indel) Y TAMBIÉN un evento
    intermedio (> max_indel, < large_del_threshold), se conservan AMBOS en
    la llamada (no se descarta el intermedio silenciosamente); la lectura
    cuenta una sola vez en el numerador de microindels.
    """
    lo = site.cut_pos - window
    hi = site.cut_pos + window
    if aln.mapq < min_mapq:
        return SiteCall(site, "NOT_EVALUABLE",
                        reason=f"MAPQ {aln.mapq} < {min_mapq}")
    span_start = aln.r_st
    span_end = m_blocks[-1][1] if m_blocks else aln.r_st
    # la ventana es CERRADA [lo, hi]: la lectura debe CONTENER la base hi;
    # r_en es exclusivo (primera base no cubierta) → exigir span_end > hi
    if not (span_start <= lo and span_end > hi):
        return SiteCall(site, "NOT_EVALUABLE",
                        reason="la lectura no atraviesa la ventana completa")

    def d_overlaps(ev):
        # ventana inclusive en el extremo derecho: una deleción que EMPIEZA
        # exactamente en hi (cut+window) aún se considera dentro
        return (ev.op == "D" and ev.ref_start <= hi
                and (ev.ref_start + ev.length) > lo)

    # Deleción estructural que afecta al sitio (conserva microindels asociados)
    large = [ev for ev in events if d_overlaps(ev)
             and ev.length >= large_del_threshold]
    if large:
        small = [ev for ev in events
                 if ev.op != "D" or not d_overlaps(ev)
                 or ev.length < large_del_threshold]
        small = [e for e in small if (e.op == "D" and d_overlaps(e))
                 or (e.op == "I" and lo <= e.ref_start <= hi)]
        return SiteCall(site, "LARGE_DEL", events=small, large_events=large,
                        reason="; ".join(f"deleción {e.length} nt" for e in large))

    # Microindels / eventos intermedios que solapan la ventana
    inside_d = [ev for ev in events if d_overlaps(ev)]
    inside_i = [ev for ev in events
                if ev.op == "I" and lo <= ev.ref_start <= hi]
    inside = inside_d + inside_i
    for ev in inside:
        if ev.op == "D" and ev.mean_q is None and fqual:
            fp = _locate_in_fragment(aln, ev.ref_start)
            if fp is not None:
                ev.mean_q = _deletion_flank_q(fqual, fp, frag_len)
    micro = [ev for ev in inside if ev.length <= max_indel]
    mid = [ev for ev in inside if max_indel < ev.length]
    if micro:
        # El microindel manda en la categoría, pero los eventos intermedios
        # (> max_indel) que solapan la ventana se CONSERVAN en la salida:
        # forman parte del alelo y no se descartan silenciosamente.
        reason = ""
        if mid:
            reason = ("además evento intermedio > %d nt: %s"
                      % (max_indel, "; ".join(
                          f"{'DEL' if e.op == 'D' else 'INS'}{e.length} nt"
                          for e in mid)))
        return SiteCall(site, "MICROINDEL", events=inside, reason=reason)
    if mid:
        return SiteCall(site, "OTRO_INDEL", events=mid,
                        reason="; ".join(
                            f"{'DEL' if e.op == 'D' else 'INS'}{e.length} nt "
                            f"(> {max_indel})" for e in mid))
    return SiteCall(site, "NO_INDEL")


# ---------------------------------------------------------------- análisis principal

def run_microindel_analysis(
    fastq_gz: Path,
    ref_fa: Path,
    guides: Dict[str, str],
    out_prefix: Path,
    window: int = 20,
    max_indel: int = 50,
    large_del_threshold: int = 100,
    min_mapq: int = 20,
    threads: int = 4,
    nuclease: str = "cas9",
    log=print,
) -> dict:
    """Ejecuta el análisis y escribe las tres tablas TSV."""
    import pandas as pd

    aligner = mp.Aligner(str(ref_fa), preset="map-ont", n_threads=threads,
                        best_n=5, bw=2000, bw_long=2000)
    if not aligner:
        raise RuntimeError(f"No se pudo indexar {ref_fa}")
    ref_name = aligner.seq_names[0]
    ref_seq = aligner.seq(ref_name).upper()

    log(f"  · localizando guías en {ref_name} ({len(ref_seq)} nt)")
    sites = locate_cut_sites(ref_seq, guides, nuclease=nuclease, log=log)
    by_id = {s.site_id: s for s in sites}

    per_read_rows = []
    allele_counter: Counter = Counter()
    allele_reads: Dict[tuple, List[str]] = defaultdict(list)
    site_stats = {s.site_id: Counter() for s in sites}
    n_reads = 0

    for name, seq, qual in mp.fastx_read(str(fastq_gz), read_comment=False):
        n_reads += 1
        hits = list(aligner.map(seq, cs=False, MD=False))
        if not hits:
            for s in sites:
                per_read_rows.append({
                    "read_id": name, "sitio": s.site_id,
                    "categoria": "NOT_EVALUABLE", "motivo": "sin alineamiento",
                    "n_eventos": 0, "eventos": "", "eventos_grandes": "",
                    "alelo": "", "max_indel_nt": 0, "calidad_media_eventos": "",
                })
                site_stats[s.site_id]["NOT_EVALUABLE"] += 1
            continue

        # Para cada sitio: el mejor alineamiento que ATRAVIESE la ventana;
        # si ninguno la atraviesa, el más largo (para justificar NOT_EVALUABLE).
        # Si VARIOS alineamientos atraviesan la ventana y discrepan sobre si
        # hay evento dentro (lecturas quiméricas con flancos duplicados),
        # el sitio se marca NOT_EVALUABLE en lugar de elegir arbitrariamente.
        for s in sites:
            lo, hi = s.cut_pos - window, s.cut_pos + window
            # ventana CERRADA [lo, hi]: r_en es exclusivo → exigir r_en > hi
            covering = [h for h in hits
                        if h.r_st <= lo and h.r_en > hi
                        and h.mapq >= min_mapq]
            discordant = False
            if len(covering) > 1:
                patrones = []
                for h in covering:
                    frag_h, fq_h = _oriented_fragment(h, seq, qual)
                    evs_h, _ = _extract_events(h, frag_h, fq_h)
                    dentro_h = [e for e in evs_h if (
                        (e.op == "D" and e.ref_start <= hi
                         and e.ref_start + e.length > lo)
                        or (e.op == "I" and lo <= e.ref_start <= hi))]
                    patrones.append(tuple(sorted(
                        (e.op, e.ref_start, e.length, e.ins_seq)
                        for e in normalize_events(dentro_h, ref_seq))))
                if len(set(patrones)) > 1:
                    discordant = True
            if discordant:
                per_read_rows.append({
                    "read_id": name, "sitio": s.site_id,
                    "categoria": "NOT_EVALUABLE",
                    "motivo": "alineamientos discordantes en la ventana: " + repr(patrones),
                    "n_eventos": 0, "eventos": "", "eventos_grandes": "",
                    "alelo": "", "max_indel_nt": 0,
                    "calidad_media_eventos": "",
                })
                site_stats[s.site_id]["NOT_EVALUABLE"] += 1
                continue
            aln = max(covering or hits, key=lambda h: h.mlen)
            frag, fqual = _oriented_fragment(aln, seq, qual)
            events, m_blocks = _extract_events(aln, frag, fqual)
            call = call_site(s, aln, events, m_blocks, fqual, len(frag),
                             window, large_del_threshold, min_mapq, max_indel)
            st = site_stats[s.site_id]
            st[call.category] += 1
            if call.category == "MICROINDEL" and call.events:
                evs = normalize_events(call.events, ref_seq)
                key = tuple((e.op, e.ref_start, e.length, e.ins_seq)
                           for e in evs)
                allele_counter[(s.site_id, key)] += 1
                allele_reads[(s.site_id, key)].append(name)
                max_len = max((e.length for e in evs), default=0)
                mqs = [e.mean_q for e in evs if e.mean_q is not None]
                per_read_rows.append({
                    "read_id": name, "sitio": s.site_id,
                    "categoria": call.category,
                    "motivo": call.reason,
                    "n_eventos": len(evs),
                    "eventos": "; ".join(str(e) for e in evs),
                    "eventos_grandes": "; ".join(str(e) for e in call.large_events),
                    "alelo": f"{s.site_id}|{key}",
                    "max_indel_nt": max_len,
                    "calidad_media_eventos": f"{sum(mqs)/len(mqs):.1f}" if mqs else "",
                })
            elif call.category in ("MICROINDEL", "LARGE_DEL", "OTRO_INDEL"):
                # LARGE_DEL / OTRO_INDEL: serializar SIEMPRE los eventos
                # (grandes y pequeños) para no perderlos silenciosamente
                evs = normalize_events(call.events, ref_seq) if call.events else []
                max_len = max((e.length for e in evs), default=0)
                mqs = [e.mean_q for e in evs if e.mean_q is not None]
                per_read_rows.append({
                    "read_id": name, "sitio": s.site_id,
                    "categoria": call.category,
                    "motivo": call.reason,
                    "n_eventos": len(evs),
                    "eventos": "; ".join(str(e) for e in evs),
                    "eventos_grandes": "; ".join(str(e) for e in call.large_events),
                    "alelo": "",
                    "max_indel_nt": max_len,
                    "calidad_media_eventos": f"{sum(mqs)/len(mqs):.1f}" if mqs else "",
                })
            else:
                per_read_rows.append({
                    "read_id": name, "sitio": s.site_id,
                    "categoria": call.category, "motivo": call.reason,
                    "n_eventos": 0, "eventos": "", "eventos_grandes": "",
                    "alelo": "", "max_indel_nt": 0, "calidad_media_eventos": "",
                })

    df_read = pd.DataFrame(per_read_rows)
    df_read.to_csv(f"{out_prefix}.microindels_per_read.tsv", sep="\t",
                   index=False)

    # --- Resumen por sitio
    rows = []
    for s in sites:
        st = site_stats[s.site_id]
        ev = st["NO_INDEL"] + st["MICROINDEL"]
        rows.append({
            "sitio": s.site_id,
            "guia": s.guide_label,
            "hebra": s.strand,
            "corte_1based": s.cut_1based,
            "lecturas_totales": n_reads,
            "sin_alineamiento_o_no_evaluables": st["NOT_EVALUABLE"],
            "sitio_eliminado_por_deleccion_grande": st["LARGE_DEL"],
            "otro_indel_intermedio": st["OTRO_INDEL"],
            "evaluables": ev,
            "sin_microindel": st["NO_INDEL"],
            "con_microindel": st["MICROINDEL"],
            "pct_microindel_sobre_evaluables": (
                f"{100 * st['MICROINDEL'] / ev:.2f}" if ev else "n/d"),
        })
    df_sum = pd.DataFrame(rows)
    df_sum.to_csv(f"{out_prefix}.microindels_summary.tsv", sep="\t",
                  index=False)

    # --- Alelos + contexto de unión RECONSTRUIDO desde la referencia
    allele_rows = []
    allele_ids: Dict[tuple, str] = {}
    for (site_id, key), n in allele_counter.most_common():
        s = by_id[site_id]
        evs = list(key)
        allele_ids[(site_id, key)] = f"{site_id}|alelo_{len(allele_ids) + 1:03d}"
        lo = max(0, s.cut_pos - 30)
        hi = min(len(ref_seq), s.cut_pos + 30)
        ref_win = ref_seq[lo:hi]
        read_win = list(ref_win)
        for op, rpos, length, ins in sorted(evs, key=lambda e: -e[1]):
            rel = rpos - lo
            if op == "D":
                a, b = max(0, rel), min(len(read_win), rel + length)
                if a < b:
                    del read_win[a:b]
            elif op == "I" and 0 <= rel <= len(read_win):
                read_win.insert(rel, f"[{ins}]")
        ev_total = site_stats[site_id]['NO_INDEL'] + site_stats[site_id]['MICROINDEL']
        allele_rows.append({
            "alelo": allele_ids[(site_id, key)],
            "sitio": site_id,
            "n_lecturas": n,
            "pct_sobre_evaluables": f"{100 * n / ev_total:.2f}" if ev_total else "n/d",
            "eventos": "; ".join(
                f"{'DEL' if op == 'D' else 'INS'}{length}nt@ref{rpos + 1}"
                for op, rpos, length, ins in evs),
            "secuencia_insercion": " | ".join(
                ins for op, rpos, length, ins in evs if op == "I"),
            "contexto_ref_(corte±30)": ref_win,
            "contexto_RECONSTRUIDO_desde_ref_(corte±30)": "".join(read_win),
            "ids_lecturas": ",".join(allele_reads[(site_id, key)][:20]),
        })
    df_alleles = pd.DataFrame(allele_rows)
    df_alleles.to_csv(f"{out_prefix}.microindels_alleles.tsv", sep="\t",
                      index=False)

    # Reetiquetar la columna "alelo" con ids cortos
    if not df_read.empty:
        df_read["alelo"] = df_read["alelo"].replace(
            {f"{sid}|{k}": aid for (sid, k), aid in allele_ids.items()})
        df_read.to_csv(f"{out_prefix}.microindels_per_read.tsv", sep="\t",
                       index=False)

    # Comprobación de coherencia numerador/denominador
    for s in sites:
        st = site_stats[s.site_id]
        total_calls = sum(st.values())
        assert total_calls == n_reads, \
            f"{s.site_id}: {total_calls} llamadas != {n_reads} lecturas"

    log(f"  · microindels: {n_reads} lecturas, {len(sites)} sitios, "
        f"{sum(1 for r in per_read_rows if r['categoria'] == 'MICROINDEL')} "
        f"llamadas con microindel")

    # Guías solicitadas que NO se localizaron en el amplicón WT: se avisa de
    # forma explícita (antes se omitían en silencio y solo quedaba un aviso
    # enterrado en el registro).
    encontradas = {s.guide_label for s in sites}
    no_localizadas = [lbl for lbl in guides if lbl not in encontradas]
    for lbl in no_localizadas:
        log(f"  [AVISO] La guía '{lbl}' NO se localiza en el amplicón WT: "
            f"se omite en el análisis de microindels. Revisa que la secuencia "
            f"corresponde al alelo salvaje (no al editado) y que no atraviesa "
            f"la unión de la deleción.")

    return {
        "n_reads": n_reads,
        "n_sites": len(sites),
        "sites": [s.site_id for s in sites],
        "guias_no_localizadas": no_localizadas,
        "summary": df_sum.to_dict("records"),
        "files": {
            "per_read": f"{out_prefix}.microindels_per_read.tsv",
            "summary": f"{out_prefix}.microindels_summary.tsv",
            "alleles": f"{out_prefix}.microindels_alleles.tsv",
        },
    }


def parse_guides(text: str) -> Dict[str, str]:
    """
    Analiza el texto del formulario/servidor: una guía por línea.
      - 'ACGT...ACGT'                  → etiqueta automática guia_1…
      - 'mi_etiqueta: ACGT...ACGT'     → etiqueta personalizada
      - 'mi_etiqueta<TAB>ACGT...ACGT'
    Acepta ARN (U) o ADN (T).
    """
    guides: Dict[str, str] = {}
    for i, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if ":" in line:
            label, seq = line.split(":", 1)
            label = label.strip() or f"guia_{i}"
        elif "\t" in line:
            label, seq = line.split("\t", 1)
            label = label.strip() or f"guia_{i}"
        else:
            label, seq = f"guia_{i}", line
        seq = _norm_guide(seq.strip().replace(" ", ""))
        if seq:
            guides[label] = seq
    return guides
