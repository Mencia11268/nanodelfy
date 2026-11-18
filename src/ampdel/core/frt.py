"""Deleción teórica y cambios locales en FRT, con exclusiones conservadoras.

Perfecta significa secuencia de unión compatible con la referencia teórica,
no identificación única del punto de corte biológico dentro de la repetición.
No infiere extremos a partir de segmentos solapados. Compatible con Windows.
"""
from __future__ import annotations
import csv
import json
from collections import Counter
from pathlib import Path
import mappy as mp
from .microindel import (
    _extract_events, _oriented_fragment, normalize_events, locate_cut_sites,
)

CLASES = (
    "DELECION_TEORICA_UNION_EXACTA",
    "DELECION_GRANDE_UNION_DIFERENTE",
    "DELECION_GRANDE_UNION_NO_RESUELTA",
    "SIN_DELECION_GRANDE", "NO_EVALUABLE",
)


def escribir_tsv(ruta, filas, campos=None):
    campos = campos or list(filas[0])
    with open(ruta, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=campos, delimiter="\t")
        w.writeheader()
        w.writerows(filas)


def leer_fasta(ruta):
    return "".join(x.strip() for x in Path(ruta).read_text().splitlines()
                   if not x.startswith(">")).upper()


def _bordes_query(hit, posicion):
    """Fronteras exactas dentro de M, nunca dentro de D ni de I."""
    r, q = hit.r_st, 0
    resultado = []
    for longitud, op in hit.cigar:
        if op in (0, 7, 8):
            if r <= posicion <= r + longitud:
                resultado.append(q + posicion - r)
            r += longitud
            q += longitud
        elif op == 1:
            q += longitud
        elif op in (2, 3):
            r += longitud
    return sorted(set(resultado))


def observacion(hit, seq, qual, ref, lo, hi):
    """Evidencia localizada [lo,hi), orientada a la referencia."""
    frag, calidad = _oriented_fragment(hit, seq, qual)
    eventos, _ = _extract_events(hit, frag, calidad)
    dentro = [e for e in eventos if
              (e.op == "D" and e.ref_start < hi and e.ref_start + e.length > lo)
              or (e.op == "I" and lo <= e.ref_start < hi)]
    normalizados = normalize_events(dentro, ref)
    patron = [(e.op, e.ref_start, e.length, e.ins_seq) for e in normalizados]
    sustituciones = []
    r, q = hit.r_st, 0
    for longitud, op in hit.cigar:
        if op in (0, 7, 8):
            for p in range(max(lo, r), min(hi, r + longitud)):
                base = frag[q + p - r].upper()
                if base != ref[p]:
                    sustituciones.append(("X", p, ref[p], base))
            r += longitud
            q += longitud
        elif op == 1:
            q += longitud
        elif op in (2, 3):
            r += longitud
    patron.extend(sustituciones)
    a, b = _bordes_query(hit, lo), _bordes_query(hit, hi)
    contexto = frag[min(a):max(b)] if a and b else ""
    return {"patron": tuple(sorted(patron)), "eventos": dentro,
            "sustituciones": sustituciones, "contexto": contexto}


def llamada_local(hits, seq, qual, ref, lo, hi, min_mapq=20):
    """Compara alelos completos, incluida la secuencia insertada."""
    elegibles = [h for h in hits if h.mapq >= min_mapq
                 and h.r_st <= lo - 10 and h.r_en >= hi + 10]
    obs = [observacion(h, seq, qual, ref, lo, hi) for h in elegibles]
    if not obs:
        return None, "sin cobertura local con ambos flancos", []
    evidencia = [list(o["patron"]) for o in obs]
    if len({o["patron"] for o in obs}) > 1:
        return None, "alelos discordantes entre alineamientos", evidencia
    i = max(range(len(elegibles)), key=lambda i: elegibles[i].mlen)
    return obs[i], "", evidencia


def _fraccion(hit, seq):
    return (hit.q_en - hit.q_st) / max(len(seq), 1)


def clasificar_estructura(seq, qual, wt, edit, hits_wt, hits_edit,
                         ds, de, min_mapq=20, umbral_grande=100):
    """Clasifica por evidencia de ambos flancos, no por longitud de lectura."""
    buenos_wt = [h for h in hits_wt if h.mapq >= min_mapq]
    buenos_ed = [h for h in hits_edit if h.mapq >= min_mapq]
    cruzan = [h for h in buenos_wt if h.r_st <= ds - 60 and h.r_en >= de + 60]
    eventos_grandes, sin_grande = [], []
    for h in cruzan:
        frag, fq = _oriented_fragment(h, seq, qual)
        eventos, _ = _extract_events(h, frag, fq)
        grandes = [e for e in eventos if e.op == "D" and e.length >= umbral_grande
                   and e.ref_start < de + 60 and e.ref_start + e.length > ds - 60]
        if grandes:
            eventos_grandes.append((h, grandes))
        else:
            sin_grande.append(h)
    jun = [h for h in buenos_ed if h.r_st <= ds - 80 and h.r_en >= ds + 80]
    mejor_wt = max((h.mlen for h in buenos_wt), default=0)
    mejor_ed = max((h.mlen for h in jun), default=0)
    evidencia_edit = bool(jun) and mejor_ed > mejor_wt * 1.05
    if sin_grande and evidencia_edit:
        return "NO_EVALUABLE", "evidencia de estructuras incompatibles", {}
    if eventos_grandes and sin_grande:
        return "NO_EVALUABLE", "retención y deleción en segmentos distintos", {}
    if eventos_grandes or evidencia_edit:
        candidatos = [h for h in jun if _fraccion(h, seq) >= .80]
        if not candidatos:
            if eventos_grandes and any(_fraccion(h, seq) >= .80
                                      for h, _ in eventos_grandes):
                return "DELECION_GRANDE_UNION_DIFERENTE", "D >=100 nt con flancos en WT", {
                    "eventos_grandes": [str(e) for _, es in eventos_grandes for e in es]}
            return "NO_EVALUABLE", "cobertura insuficiente o estructura compleja", {}
        o, motivo, evidencia = llamada_local(
            candidatos, seq, qual, edit, ds - 30, ds + 30, min_mapq)
        if o is None:
            if "discordantes" in motivo:
                return "NO_EVALUABLE", motivo, {"patrones_union": evidencia}
            return "DELECION_GRANDE_UNION_NO_RESUELTA", motivo, {}
        info = {
            "contexto_union_observado": o["contexto"],
            "contexto_union_teorico": edit[ds - 30:ds + 30],
            "patron_union": list(o["patron"]),
            "nota": "compatible; extremos no identificables de forma única en FRT",
        }
        if not o["patron"] and o["contexto"] == edit[ds - 30:ds + 30]:
            return "DELECION_TEORICA_UNION_EXACTA", "unión exacta en 60 nt", info
        if not o["contexto"] and not o["patron"]:
            return "DELECION_GRANDE_UNION_NO_RESUELTA", "borde local dentro de hueco", info
        # Una D que cruza el borde de la ventana impide extraer 60 nt,
        # pero demuestra que la unión no es idéntica a la teórica.
        # No inventamos coordenadas de unión ni la llamamos perfecta.
        return "DELECION_GRANDE_UNION_DIFERENTE", "diferencias en unión de 60 nt", info
    if sin_grande:
        if not any(_fraccion(h, seq) >= .80 for h in sin_grande):
            return "NO_EVALUABLE", "posible concatémero o lectura compleja", {}
        return "SIN_DELECION_GRANDE", "retención con ambos FRT flanqueados", {}
    return "NO_EVALUABLE", "no distingue retención del intervalo y deleción", {}


def _categoria_local(o):
    if o is None:
        return "NO_EVALUABLE"
    ev, sub = o["eventos"], o["sustituciones"]
    if any(e.op == "D" and e.length >= 100 for e in ev):
        return "DELECION_LOCAL_GRANDE"
    if any(e.length > 50 for e in ev):
        return "INDEL_INTERMEDIO"
    if ev and sub:
        return "INDEL_Y_SUSTITUCIONES"
    if ev:
        return "INDEL_PEQUENO"
    if sub:
        return "SOLO_SUSTITUCIONES"
    return "SIN_CAMBIOS"


def analizar_frt(fastq, wt_fa, edit_fa, prefijo, guia="TCCTATTCTCTAGAAAGTAT",
                 min_mapq=20, ds=None, de=None):
    wt, edit = leer_fasta(wt_fa), leer_fasta(edit_fa)
    sitios = locate_cut_sites(wt, {"FRT": guia}, log=lambda _: None)
    assert len(sitios) == 2 and all(s.pam_found for s in sitios)
    ds = sitios[0].cut_pos if ds is None else ds
    de = sitios[1].cut_pos if de is None else de
    assert wt[:ds] + wt[de:] == edit, "Referencia teórica incompatible con esos cortes"
    intervalos = [(s.protospacer_start - 4, s.protospacer_start + 30) for s in sitios]
    assert wt[slice(*intervalos[0])] == wt[slice(*intervalos[1])]
    aw = mp.Aligner(str(wt_fa), preset="map-ont", best_n=5, bw=2000, bw_long=2000)
    ae = mp.Aligner(str(edit_fa), preset="map-ont", best_n=5, bw=2000, bw_long=2000)
    filas, locales = [], []
    for rid, seq, qual in mp.fastx_read(str(fastq)):
        hw, he = list(aw.map(seq)), list(ae.map(seq))
        clase, motivo, info = clasificar_estructura(
            seq, qual, wt, edit, hw, he, ds, de, min_mapq)
        fila = {"read_id": rid, "longitud": len(seq), "categoria": clase,
                "motivo": motivo, "evidencia": json.dumps(info, ensure_ascii=False),
                "FRT1": "NO_EVALUABLE", "FRT2": "NO_EVALUABLE"}
        for i, (lo, hi) in enumerate(intervalos, 1):
            if clase == "SIN_DELECION_GRANDE":
                o, razon, evid = llamada_local(hw, seq, qual, wt, lo, hi, min_mapq)
                cat = _categoria_local(o)
            else:
                o, razon, evid = None, "fuera del grupo sin deleción grande", []
                cat = "NO_EVALUABLE"
            fila[f"FRT{i}"] = cat
            fila[f"microindel_FRT{i}"] = bool(o and any(e.length <= 50 for e in o["eventos"]))
            fila[f"sustituciones_FRT{i}"] = bool(o and o["sustituciones"])
            locales.append({
                "read_id": rid, "sitio": f"FRT{i}", "inicio_1based": lo + 1,
                "fin_1based": hi, "categoria_estructural": clase,
                "categoria_local": cat, "motivo": razon,
                "contexto_observado": o["contexto"] if o else "",
                "contexto_referencia": wt[lo:hi],
                "eventos": "; ".join(str(e) for e in o["eventos"]) if o else "",
                "tiene_microindel_1_50nt": fila[f"microindel_FRT{i}"],
                "tiene_sustituciones": fila[f"sustituciones_FRT{i}"],
                "sustituciones": json.dumps(o["sustituciones"] if o else []),
                "patrones_alineamientos": json.dumps(evid),
            })
        filas.append(fila)
    conteo = Counter(f["categoria"] for f in filas)
    total = len(filas)
    informativas = total - conteo["NO_EVALUABLE"]
    grandes = sum(conteo[c] for c in CLASES if c.startswith("DELECION_"))
    resumen = []
    for cat in CLASES + ("TOTAL_DELECIONES_GRANDES", "INFORMATIVAS", "FILTRADAS"):
        n = (grandes if cat == "TOTAL_DELECIONES_GRANDES" else
             informativas if cat == "INFORMATIVAS" else
             total if cat == "FILTRADAS" else conteo[cat])
        resumen.append({
            "categoria": cat, "lecturas": n, "denominador_filtradas": total,
            "pct_filtradas": round(100*n/total, 3) if total else "",
            "denominador_informativas": informativas,
            "pct_informativas": (round(100*n/informativas, 3)
                                if informativas and cat not in ("NO_EVALUABLE", "FILTRADAS") else ""),
        })
    resumen_local = []
    for sitio in ("FRT1", "FRT2", "AMBOS_FRT"):
        grupo = [f for f in filas if f["categoria"] == "SIN_DELECION_GRANDE"]
        if sitio == "AMBOS_FRT":
            ev = [f for f in grupo if all(f[s] != "NO_EVALUABLE" for s in ("FRT1", "FRT2"))]
            cambios = [f for f in ev if any(f[s] != "SIN_CAMBIOS" for s in ("FRT1", "FRT2"))]
            indels = [f for f in ev if any(f[f"microindel_{s}"] for s in ("FRT1", "FRT2"))]
            sub = [f for f in ev if any(f[f"sustituciones_{s}"] for s in ("FRT1", "FRT2"))]
        else:
            ev = [f for f in grupo if f[sitio] != "NO_EVALUABLE"]
            cambios = [f for f in ev if f[sitio] != "SIN_CAMBIOS"]
            indels = [f for f in ev if f[f"microindel_{sitio}"]]
            sub = [f for f in ev if f[f"sustituciones_{sitio}"]]
        den = len(ev)
        resumen_local.append({
            "sitio": sitio, "sin_delecion_grande": len(grupo),
            "evaluables": den, "con_cambios": len(cambios),
            "con_indels_1_50nt": len(indels), "con_sustituciones": len(sub),
            "pct_cambios": round(100*len(cambios)/den, 3) if den else "",
            "pct_indels": round(100*len(indels)/den, 3) if den else "",
            "pct_sustituciones": round(100*len(sub)/den, 3) if den else "",
        })
    prefijo = str(prefijo)
    escribir_tsv(prefijo + ".frt_estructura_lecturas.tsv", filas)
    escribir_tsv(prefijo + ".frt_estructura_resumen.tsv", resumen)
    escribir_tsv(prefijo + ".frt_sin_delecion_lecturas.tsv", locales)
    escribir_tsv(prefijo + ".frt_sin_delecion_resumen.tsv", resumen_local)
    return {"estructura": resumen, "local": resumen_local, "n": total}
