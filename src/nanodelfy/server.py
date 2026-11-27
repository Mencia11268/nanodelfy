"""Servidor FastAPI: rutas web + API JSON para lanzar y consultar trabajos."""
from __future__ import annotations
import os, re, uuid, shutil, traceback, threading, json
from pathlib import Path
from datetime import datetime
from typing import Dict, Any

from fastapi import FastAPI, Request, UploadFile, File, Form, HTTPException
from fastapi.responses import HTMLResponse, JSONResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from .core.pipeline import run_full_pipeline


# --------- rutas de almacenamiento -----------------------------------------
BASE_DIR = Path(__file__).resolve().parent
WORK_DIR = Path(os.environ.get("NANODEL_WORK_DIR",
                               Path.cwd() / "uploads")).resolve()
WORK_DIR.mkdir(parents=True, exist_ok=True)

app = FastAPI(title="nanodelfy")
templates = Jinja2Templates(directory=str(BASE_DIR / "templates"))
app.mount("/static", StaticFiles(directory=str(BASE_DIR / "static")), name="static")

# Estado en memoria de los trabajos lanzados (no persistente)
JOBS: Dict[str, Dict[str, Any]] = {}


def _safe_sample_name(name: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    return name or "sample"


# ============================================================
# Rutas web
# ============================================================

@app.get("/", response_class=HTMLResponse)
def index(request: Request):
    return templates.TemplateResponse(request, "index.html", {})


@app.get("/job/{job_id}", response_class=HTMLResponse)
def job_page(request: Request, job_id: str):
    if job_id not in JOBS:
        raise HTTPException(404, "Trabajo no encontrado")
    return templates.TemplateResponse(
        request, "job.html",
        {"job_id": job_id, "job": JOBS[job_id]},
    )


# ============================================================
# API
# ============================================================

@app.post("/api/run")
async def api_run(
    fastq: list[UploadFile] = File(...),
    wt_fasta: UploadFile = File(...),
    sample_name: str = Form(""),
    deletion_start: str = Form(""),
    deletion_end: str = Form(""),
    min_quality: int = Form(10),
    min_length: int = Form(800),
    max_length: int = Form(4000),
    breakpoint_tolerance: int = Form(20),
    min_sv_deletion: int = Form(1500),
    min_edit_ref_coverage: float = Form(0.90),
    min_mapq: int = Form(20),
    wt_min: int = Form(2800),
    wt_max: int = Form(3500),
    nhej_search_window: int = Form(500),
    min_junction_insertion: int = Form(5),
    guides: str = Form(""),
    nuclease: str = Form("cas9"),
    guide_max_mismatch: int = Form(1),
    snp_min_alt_freq: float = Form(0.20),
    snp_min_depth: int = Form(10),
    calibration: bool = Form(True),
    microindel_window: int = Form(20),
    microindel_max_indel: int = Form(50),
    microindel_large_del_threshold: int = Form(100),
):
    # Validaciones rápidas
    def _coord(valor: str, nombre: str):
        valor = (valor or "").strip()
        if not valor:
            return None
        try:
            n = int(valor)
        except ValueError:
            raise HTTPException(400, f"{nombre} debe ser un número entero.")
        if n < 1:
            raise HTTPException(400, f"{nombre} debe ser >= 1 (1-based).")
        return n

    del_start = _coord(deletion_start, "Inicio de la deleción")
    del_end = _coord(deletion_end, "Fin de la deleción")
    if (del_start is None) != (del_end is None):
        raise HTTPException(400, "Indica ambas coordenadas de la deleción "
                                  "o ninguna (se calculan desde las guías).")
    if del_start is not None and del_end <= del_start:
        raise HTTPException(400, "Coordenadas inválidas: end debe ser mayor "
                                 "que start (1-based).")
    if del_start is None and not (guides or "").strip():
        raise HTTPException(400, "Sin coordenadas de la deleción necesitas "
                                  "al menos dos guías para calcularlas "
                                  "automáticamente.")

    job_id = datetime.now().strftime("%Y%m%d_%H%M%S_") + uuid.uuid4().hex[:6]
    job_dir = WORK_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    def _nombre_desde_archivo(filename: str) -> str:
        """Nombre de muestra a partir del archivo (sin extensión fastq)."""
        nombre = (filename or "").rsplit("/", 1)[-1].rsplit("\\", 1)[-1]
        for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if nombre.lower().endswith(suf):
                nombre = nombre[: -len(suf)]
                break
        return _safe_sample_name(nombre)

    # Varios FASTQ en la misma petición: una muestra por archivo, con las
    # mismas condiciones. El nombre de cada muestra sale de su archivo.
    muestras = []          # [(nombre, ruta_fastq)]
    vistos = {}
    for f in fastq:
        fastq_name = f.filename or "input.fastq.gz"
        if not fastq_name.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
            raise HTTPException(400, "El archivo de lecturas debe ser FASTQ "
                                      "(.fastq/.fastq.gz/.fq/.fq.gz).")
        nombre = _nombre_desde_archivo(fastq_name)
        vistos[nombre] = vistos.get(nombre, 0) + 1
        if vistos[nombre] > 1:
            nombre = f"{nombre}_{vistos[nombre]}"
        fastq_path = job_dir / f"input_{len(muestras)}_{nombre}.fastq"
        with open(fastq_path, "wb") as out:
            shutil.copyfileobj(f.file, out)
        muestras.append((nombre, fastq_path))
    if not muestras:
        raise HTTPException(400, "Sube al menos un archivo FASTQ.")

    wt_path = job_dir / "amplicon_WT_input.fa"
    with open(wt_path, "wb") as f:
        shutil.copyfileobj(wt_fasta.file, f)

    # Etiqueta del análisis (cabecera): la escrita en el formulario, o el
    # nombre de la única muestra / "varias muestras" si se dejó vacía
    label = _safe_sample_name(sample_name) if (sample_name or "").strip() else (
        muestras[0][0] if len(muestras) == 1 else "varias_muestras")
    sample_name = label

    params = {
        "deletion_start": del_start,
        "deletion_end": del_end,
        "min_quality": min_quality,
        "min_length": min_length,
        "max_length": max_length,
        "breakpoint_tolerance": breakpoint_tolerance,
        "min_sv_deletion": min_sv_deletion,
        "min_edit_ref_coverage": min_edit_ref_coverage,
        "min_mapq": min_mapq,
        "wt_min": wt_min,
        "wt_max": wt_max,
        "nhej_search_window": nhej_search_window,
        "min_junction_insertion": min_junction_insertion,
        "guides": guides,
        "nuclease": nuclease,
        "guide_max_mismatch": guide_max_mismatch,
        "snp_min_alt_freq": snp_min_alt_freq,
        "snp_min_depth": snp_min_depth,
        "calibration": calibration,
        "microindel_window": microindel_window,
        "microindel_max_indel": microindel_max_indel,
        "microindel_large_del_threshold": microindel_large_del_threshold,
    }

    job = {
        "id": job_id,
        "sample_name": sample_name,
        "status": "running",
        "progress": 0,
        "stage": "Iniciando",
        "log": [],
        "result": None,
        "error": None,
        "params": params,
        "job_dir": str(job_dir),
        "created_at": datetime.now().isoformat(timespec="seconds"),
    }
    JOBS[job_id] = job

    def _log(msg: str):
        job["log"].append(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")

    def _set_progress(p: int, s: str):
        job["progress"] = p
        job["stage"] = s

    def _worker():
        n = len(muestras)
        try:
            resultados = []
            for i, (nombre, fastq_path) in enumerate(muestras):
                def _log_m(msg, _n=nombre):
                    _log(f"[{_n}] {msg}")

                def _set_progress_m(p, s, _i=i, _n=nombre, _n_total=n):
                    # progreso global: fracción de muestras terminadas +
                    # fracción de la muestra en curso
                    job["progress"] = int(100 * (_i + p / 100) / _n_total)
                    job["stage"] = f"Muestra {_i + 1}/{_n_total} · {_n}: {s}"

                _log(f"=== Muestra {i + 1}/{n}: {nombre} ===")
                sub = run_full_pipeline(
                    job_dir=job_dir / "samples" / f"{i:02d}_{nombre}",
                    fastq_path=fastq_path,
                    wt_fasta_path=wt_path,
                    sample_name=nombre,
                    params=params,
                    log=_log_m,
                    set_progress=_set_progress_m,
                )
                resultados.append({"name": nombre, "result": sub})
            job["result"] = {"samples": resultados}
            job["status"] = "done"
            job["progress"] = 100
            job["stage"] = "Completado"
        except Exception as e:
            job["status"] = "error"
            job["error"] = f"{type(e).__name__}: {e}"
            job["log"].append("ERROR: " + job["error"])
            job["log"].append(traceback.format_exc())

    threading.Thread(target=_worker, daemon=True).start()
    return {"job_id": job_id}


@app.get("/api/job/{job_id}")
def api_job(job_id: str):
    if job_id not in JOBS:
        raise HTTPException(404, "Trabajo no encontrado")
    return JOBS[job_id]


@app.get("/api/download/{job_id}/{kind}")
def api_download(job_id: str, kind: str):
    job = JOBS.get(job_id)
    if not job:
        raise HTTPException(404, "Trabajo no encontrado")
    if job["status"] != "done":
        raise HTTPException(400, "El trabajo aún no ha terminado")

    if kind == "all":
        """Zip con todos los resultados del trabajo (una sola descarga)."""
        import zipfile
        job_dir = Path(job["job_dir"])
        zip_path = job_dir / f"{job['sample_name']}_resultados.zip"
        if not zip_path.exists():
            # nuevo formato: una subcarpeta por muestra bajo samples/;
            # formato antiguo: todo en results/
            muestras_dir = job_dir / "samples"
            with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
                raiz = muestras_dir if muestras_dir.is_dir() else job_dir / "results"
                for p in sorted(raiz.rglob("*")):
                    if p.is_file():
                        zf.write(p, p.relative_to(raiz.parent))
        return FileResponse(zip_path, filename=zip_path.name,
                            media_type="application/zip")

    # nuevo formato multi-muestra: kind = "s<índice>_<tipo>"
    result = job["result"] or {}
    m = re.match(r"^s(\d+)_(.+)$", kind)
    if m and isinstance(result.get("samples"), list):
        idx, sub_kind = int(m.group(1)), m.group(2)
        try:
            files = result["samples"][idx]["result"]["files"]
        except (IndexError, KeyError, TypeError):
            raise HTTPException(404, f"Muestra desconocida en el trabajo: {kind}")
    elif isinstance(result.get("samples"), list):
        # formato nuevo sin índice: no válido
        raise HTTPException(404, f"Tipo de archivo desconocido: {kind} "
                                 "(indica la muestra: s0_, s1_, ...)")
    else:
        files = result.get("files", {})
        sub_kind = kind
    if sub_kind not in files:
        raise HTTPException(404, f"Tipo de archivo desconocido: {kind}")
    path = Path(files[sub_kind])
    if not path.exists():
        raise HTTPException(404, "Archivo no encontrado en disco")
    media = "application/octet-stream"
    if path.suffix == ".tsv":
        media = "text/tab-separated-values"
    elif path.suffix == ".xlsx":
        media = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    elif path.suffix == ".zip":
        media = "application/zip"
    return FileResponse(path, filename=path.name, media_type=media)
