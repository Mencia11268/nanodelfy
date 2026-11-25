"""Orquestador: lee inputs, ejecuta filtrado + alineamiento + clasificación."""
from __future__ import annotations
import shutil, zipfile, traceback
from pathlib import Path
from typing import Callable

from .refs import build_edited_reference
from .filter_align import filter_reads, align_minimap2
from .classify import Params, run_classification
from .microindel import run_microindel_analysis, parse_guides
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

    # ------- 1) Construir referencia editada
    set_progress(10, "Construyendo referencia editada")
    log("[1/4] Construyendo referencia editada...")
    wt_fa = job_dir / "amplicon_WT.fa"
    edit_fa = job_dir / "amplicon_EDIT.fa"
    shutil.copy(wt_fasta_path, wt_fa)
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
    # Parche R3EB_B: secuencia WT para reconstruir uniones desde soft-clips
    from .refs import parse_first_fasta_record
    wt_seq = str(parse_first_fasta_record(wt_fa).seq).upper()
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
    guides_text = params.get("guides", "")
    if guides_text and guides_text.strip():
        set_progress(90, "Analizando microindels en cortes de guía")
        log("[5/5] Microindels alrededor de los cortes de guía...")
        guides = parse_guides(guides_text)
        if guides:
            micro_stats = run_microindel_analysis(
                filt_fq, wt_fa, guides,
                out_prefix=out_dir / sample_name,
                window=params.get("microindel_window", 20),
                max_indel=params.get("microindel_max_indel", 50),
                large_del_threshold=params.get("microindel_large_del_threshold", 100),
                min_mapq=params.get("min_mapq", 20),
                nuclease=params.get("nuclease", "cas9"),
                guide_max_mismatch=params.get("guide_max_mismatch", 1),
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
        })
    if frt_stats:
        files.update({
            "frt_estructura": str(out_dir / f"{sample_name}.frt_estructura_resumen.tsv"),
            "frt_por_lectura": str(out_dir / f"{sample_name}.frt_estructura_lecturas.tsv"),
            "frt_sin_delecion": str(out_dir / f"{sample_name}.frt_sin_delecion_resumen.tsv"),
            "frt_local_por_lectura": str(out_dir / f"{sample_name}.frt_sin_delecion_lecturas.tsv"),
        })
    return {
        "ref_stats": ref_stats,
        "filt_stats": filt_stats,
        "classification": cls_stats,
        "microindels": micro_stats,
        "frt": frt_stats,
        "files": files,
    }
