"""Orquestador: lee inputs, ejecuta filtrado + alineamiento + clasificación."""
from __future__ import annotations
import shutil, zipfile, traceback
from pathlib import Path
from typing import Callable

from .refs import build_edited_reference, parse_first_fasta_record
from .filter_align import filter_reads, align_minimap2
from .classify import Params, run_classification
from .microindel import (run_microindel_analysis, parse_guides,
                        locate_cut_sites, _allele_segments)
from .frt import analizar_frt


def run_full_pipeline(
    job_dir: Path,
    fastq_path: Path,
    wt_fasta_path: Path,
    sample_name: str,
    params: dict,
    log: Callable[[str], None] = print,
    set_progress: Callable[[int, str], None] = lambda p, s: None,
) -> dict:
    """
    Ejecuta el pipeline en job_dir y devuelve un dict con rutas de salida y stats.
    """
    job_dir = Path(job_dir)
    out_dir = job_dir / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------- 0) Coordenadas de la deleción: explícitas o desde las guías
    wt_fa = job_dir / "amplicon_WT.fa"
    edit_fa = job_dir / "amplicon_EDIT.fa"
    shutil.copy(wt_fasta_path, wt_fa)
    wt_seq = str(parse_first_fasta_record(wt_fa).seq).upper()

    if not params.get("deletion_start") or not params.get("deletion_end"):
        guides_text = params.get("guides", "") or ""
        guides = parse_guides(guides_text) if guides_text.strip() else {}
        if len(guides) < 2:
            raise ValueError(
                "Faltan las coordenadas de la deleción (deletion_start/"
                "deletion_end). Indícalas, o introduce al menos dos guías "
                "para calcularlas de sus cortes.")
        log("Coordenadas de la deleción no indicadas: se calculan de los "
            "cortes de las guías...")
        sitios = locate_cut_sites(
            wt_seq, guides,
            nuclease=params.get("nuclease", "cas9"),
            max_mismatch=params.get("guide_max_mismatch", 1), log=log)
        pares = [(a, b) for a in sitios for b in sitios
                 if a.guide_label != b.guide_label and a.cut_pos < b.cut_pos]
        if not pares:
            raise ValueError(
                "No se localizan cortes de dos guías distintas en el "
                "amplicón WT; no se puede derivar la deleción teórica.")
        lo, hi = max(pares, key=lambda p: p[1].cut_pos - p[0].cut_pos)
        # cut_pos es 0-based "entre bases": la deleción abarca las bases
        # [lo.cut_pos, hi.cut_pos) → 1-based inclusiva lo.cut_pos+1..hi.cut_pos
        params = {**params,
                  "deletion_start": lo.cut_pos + 1,
                  "deletion_end": hi.cut_pos}
        log(f"  Deleción derivada de las guías ({lo.nuclease}): "
            f"{hi.cut_pos - lo.cut_pos} nt, coordenadas "
            f"{params['deletion_start']}-{params['deletion_end']} (1-based)")

    # ------- 1) Construir referencia editada
    set_progress(10, "Construyendo referencia editada")
    log("[1/4] Construyendo referencia editada...")
    ref_stats = build_edited_reference(
        wt_fa, edit_fa,
        deletion_start=params["deletion_start"],
        deletion_end=params["deletion_end"],
    )
    log(f"  WT: {ref_stats['wt_length']} nt · EDIT: {ref_stats['edit_length']} nt "
        f"· deleción: {ref_stats['deletion_size']} nt")

    # ------- 2) Filtrado
    set_progress(25, "Filtrando lecturas")
    log("[2/4] Filtrando por calidad y longitud...")
    filt_fq = out_dir / f"{sample_name}.filt.fastq.gz"
    filt_stats = filter_reads(
        fastq_path, filt_fq,
        min_quality=params.get("min_quality", 10),
        min_length=params.get("min_length", 200),
        max_length=params.get("max_length", 5000),
        log=log,
    )

    # ------- 3) Alineamiento doble
    set_progress(45, "Alineando contra REF_WT")
    log("[3/4] Alineando contra ambas referencias...")
    bam_wt = out_dir / f"{sample_name}.WT.bam"
    bam_edit = out_dir / f"{sample_name}.EDIT.bam"
    wt_records = align_minimap2(filt_fq, wt_fa, bam_wt, log=log)
    set_progress(65, "Alineando contra REF_EDIT")
    edit_records = align_minimap2(filt_fq, edit_fa, bam_edit, log=log)

    # ------- 4) Clasificación
    set_progress(85, "Clasificando lecturas")
    log("[4/4] Clasificando lecturas...")
    # (wt_seq se leyó en el paso 0; parche R3EB_B: reconstruir uniones
    # desde soft-clips con la secuencia WT)
    p = Params(
        deletion_start=params["deletion_start"],
        deletion_end=params["deletion_end"],
        wt_length=ref_stats["wt_length"],
        wt_seq=wt_seq,
        breakpoint_tolerance=params.get("breakpoint_tolerance", 50),
        min_sv_deletion=params.get("min_sv_deletion", 1500),
        min_edit_ref_coverage=params.get("min_edit_ref_coverage", 0.80),
        min_mapq=params.get("min_mapq", 20),
        wt_min=params.get("wt_min"),
        wt_max=params.get("wt_max"),
        nhej_search_window=params.get("nhej_search_window", 500),
        min_junction_insertion=params.get("min_junction_insertion", 5),
        calibration=params.get("calibration", True),
    )
    cls_stats = run_classification(
        wt_records, edit_records, p,
        out_tsv=out_dir / f"{sample_name}.read_classification.tsv",
        out_summary=out_dir / f"{sample_name}.summary.tsv",
        out_xlsx=out_dir / f"{sample_name}.resultados.xlsx",
    )

    # ------- 5) Microindels alrededor de los cortes de guía (opcional)
    micro_stats = None
    teorica = None
    guides_text = params.get("guides", "")
    if guides_text and guides_text.strip():
        set_progress(90, "Analizando microindels en cortes de guía")
        log("[5/5] Microindels alrededor de los cortes de guía...")
        guides = parse_guides(guides_text)
        if guides:
            # Deleción teórica entre los cortes de dos guías distintas,
            # con la regla de corte de la nucleasa elegida (cas9: corte
            # romo; cas12a: centro de la región de corte escalonado)
            try:
                sitios = locate_cut_sites(
                    wt_seq, guides,
                    nuclease=params.get("nuclease", "cas9"),
                    max_mismatch=params.get("guide_max_mismatch", 1),
                    log=log)
                pares = [(a, b) for a in sitios for b in sitios
                         if a.guide_label != b.guide_label
                         and a.cut_pos < b.cut_pos]
                if pares:
                    lo, hi = max(pares,
                                key=lambda p: p[1].cut_pos - p[0].cut_pos)
                    teorica = {
                        "start": lo.cut_pos, "end": hi.cut_pos,
                        "size": hi.cut_pos - lo.cut_pos,
                        "guia5": lo.guide_label, "guia3": hi.guide_label,
                        "nuclease": lo.nuclease,
                        "corte5": lo.cut_1based, "corte3": hi.cut_1based,
                    }
                    log(f"  Deleción teórica ({lo.nuclease}): "
                        f"{hi.cut_pos - lo.cut_pos} nt entre {lo.site_id} "
                        f"y {hi.site_id}")
            except ValueError as e:
                log(f"  [AVISO] deleción teórica no calculada: {e}")
            micro_stats = run_microindel_analysis(
                filt_fq, wt_fa, guides,
                out_prefix=out_dir / sample_name,
                window=params.get("microindel_window", 20),
                max_indel=params.get("microindel_max_indel", 50),
                large_del_threshold=params.get("microindel_large_del_threshold", 100),
                min_mapq=params.get("min_mapq", 20),
                nuclease=params.get("nuclease", "cas9"),
                guide_max_mismatch=params.get("guide_max_mismatch", 1),
                snp_min_alt_freq=params.get("snp_min_alt_freq", 0.20),
                snp_min_depth=params.get("snp_min_depth", 10),
                log=log,
            )
        else:
            log("  [AVISO] No se pudo interpretar ninguna guía; se omite el paso")

    # Módulo FRT opcional (herencia del análisis MmC7/FRT): análisis conservador
    # de la región FRT, distinto del clasificador estructural histórico.
    # Activación:
    #   - params["frt_analysis"] es True/False: explícito.
    #   - si no está definido ("auto"): se activa si alguna guía es la FRT canónica
    #     (comportamiento histórico, para los presets FRT existentes).
    # La guía del módulo se puede cambiar con params["frt_guide"].
    frt_stats = None
    _frt_mode = params.get("frt_analysis", None)
    _frt_guide = params.get("frt_guide", "TCCTATTCTCTAGAAAGTAT")
    _activate_frt = bool(_frt_mode) if _frt_mode is not None else (
        guides_text and any(g.replace("U", "T").upper() == _frt_guide
                            for g in parse_guides(guides_text).values()))
    if _activate_frt:
        try:
            frt_stats = analizar_frt(
                filt_fq, wt_fa, edit_fa, out_dir / sample_name,
                guia=_frt_guide,
                ds=params["deletion_start"] - 1, de=params["deletion_end"],
                min_mapq=params["min_mapq"])
        except AssertionError as e:
            log(f"  [AVISO] referencias no compatibles con el módulo específico FRT: {e}")

    # ------- ZIP para IGV
    set_progress(95, "Empaquetando BAMs para IGV")
    log("Empaquetando BAMs para IGV...")
    igv_zip = out_dir / f"{sample_name}.BAMs_para_IGV.zip"
    with zipfile.ZipFile(igv_zip, "w", zipfile.ZIP_STORED) as zf:
        for f in [wt_fa, wt_fa.with_suffix(".fa.fai"),
                  edit_fa, edit_fa.with_suffix(".fa.fai"),
                  bam_wt, Path(str(bam_wt) + ".bai"),
                  bam_edit, Path(str(bam_edit) + ".bai")]:
            if f.exists():
                zf.write(f, arcname=f.name)

    set_progress(100, "Listo")
    log("Pipeline completado.")
    files = {
        "summary_tsv":  str(out_dir / f"{sample_name}.summary.tsv"),
        "perread_tsv":  str(out_dir / f"{sample_name}.read_classification.tsv"),
        "excel":        str(out_dir / f"{sample_name}.resultados.xlsx"),
        "igv_zip":      str(igv_zip),
        "bam_wt":       str(bam_wt),
        "bam_edit":     str(bam_edit),
        "ref_wt":       str(wt_fa),
        "ref_edit":     str(edit_fa),
    }
    if micro_stats:
        files.update({
            "micro_per_read": micro_stats["files"]["per_read"],
            "micro_summary":  micro_stats["files"]["summary"],
            "micro_alleles":  micro_stats["files"]["alleles"],
            **({"micro_haplotipos": micro_stats["files"]["haplotipos"]}
               if micro_stats["files"].get("haplotipos") else {}),
            **({"micro_fondos": micro_stats["files"]["fondos"]}
               if micro_stats["files"].get("fondos") else {}),
        })
    if frt_stats:
        files.update({
            "frt_estructura": str(out_dir / f"{sample_name}.frt_estructura_resumen.tsv"),
            "frt_por_lectura": str(out_dir / f"{sample_name}.frt_estructura_lecturas.tsv"),
            "frt_sin_delecion": str(out_dir / f"{sample_name}.frt_sin_delecion_resumen.tsv"),
            "frt_local_por_lectura": str(out_dir / f"{sample_name}.frt_sin_delecion_lecturas.tsv"),
        })
    # --- Mapa de deleciones recurrentes (para el gráfico de lecturas web):
    # coordenadas 0-based semiabiertas, cortes de guía y deleción esperada.
    deletion_view = None
    if cls_stats and cls_stats.get("top_deletions"):
        cortes = []
        if micro_stats:
            for r in micro_stats.get("summary", []):
                cp = r.get("corte_pos_0based")
                if isinstance(cp, int):
                    cortes.append({"sitio": r["sitio"], "pos": cp})

        def _guia_cercana(pos: int):
            mejor, d = None, 60          # ±60 nt alrededor del corte
            for c in cortes:
                dd = abs(c["pos"] - pos)
                if dd <= d:
                    mejor, d = c["sitio"], dd
            return mejor

        total_reads = cls_stats.get("total_reads") or 1
        _FLANK = 30                 # nt de contexto a cada lado de la unión

        def _flancos(d: dict) -> dict:
            """Segmentos de secuencia observada alrededor de la unión.

            Usa el patrón de unión observado (cluster de la D en el CIGAR):
            intervalo [c_start, c_end) para la barra, inserción de la unión
            (bases de lectura dentro del cluster, azul) y eventos D/I de los
            flancos (tachadas en rojo / azul). Sin patrón (deleción inferida
            sin D completa en el alineamiento) se reconstruye desde la
            referencia.
            """
            pat = d.get("pattern")
            if pat and len(pat) == 4:
                c_start, c_end, ins_union, flanks = pat
                c_start, c_end = int(c_start), int(c_end)
                s0 = max(0, c_start - _FLANK)
                e1 = min(len(wt_seq), c_end + _FLANK)
                evs_f5, evs_f3 = [], []
                for op, rel, _ln, ins in flanks:
                    rpos = c_start + rel
                    ev = (op, rpos, _ln, ins or "")
                    if rpos < c_start:
                        evs_f5.append(ev)
                    elif rpos >= c_end:
                        evs_f3.append(ev)
                f5 = _allele_segments(wt_seq, s0, c_start, evs_f5) or \
                    [{"t": "m", "seq": wt_seq[s0:c_start]}]
                f3 = _allele_segments(wt_seq, c_end, e1, evs_f3) or \
                    [{"t": "m", "seq": wt_seq[c_end:e1]}]
                # el intervalo del cluster (barra + flancos) manda sobre la
                # media de breakpoints al mostrar la unión
                return {"start": c_start, "end": c_end,
                        "length": c_end - c_start,
                        "flanco5": f5, "flanco3": f3,
                        "ins_union": ins_union or ""}
            start, end = d["start"], d["end"]
            s0 = max(0, start - _FLANK)
            e1 = min(len(wt_seq), end + _FLANK)
            return {"flanco5": [{"t": "m", "seq": wt_seq[s0:start]}],
                    "flanco3": [{"t": "m", "seq": wt_seq[end:e1]}],
                    "ins_union": ""}

        deleciones = []
        for d in cls_stats["top_deletions"]:
            # sin unión observada (deleción no completa en el alineamiento,
            # p. ej. lecturas parciales) no se puede dibujar la variante:
            # fuera del gráfico (siguen contando en la clasificación)
            if not d.get("pattern"):
                continue
            item = {**d, **_flancos(d)}      # puede ajustar start/end/length
            g5, g3 = _guia_cercana(item["start"]), _guia_cercana(item["end"])
            item["pct"] = round(100.0 * d["n"] / total_reads, 2)
            item["desde_guia"] = g5 or ""
            item["hasta_guia"] = g3 or ""
            item["etiqueta"] = " → ".join(x for x in (g5, g3) if x) or \
                f"del {item['start'] + 1}-{item['end']}"
            deleciones.append(item)
            for k in ("pattern", "pattern_key"):   # internos: no van al JSON
                deleciones[-1].pop(k, None)
        deletion_view = {
            "ref_len": ref_stats["wt_length"] if ref_stats else None,
            "esperada": {"start": params["deletion_start"] - 1,
                        "end": params["deletion_end"]},
            "teorica": teorica,
            "cortes": cortes,
            "deleciones": deleciones,
        }

    return {
        "ref_stats": ref_stats,
        "filt_stats": filt_stats,
        "classification": cls_stats,
        "microindels": micro_stats,
        "frt": frt_stats,
        "deletion_view": deletion_view,
        "files": files,
    }
