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
    fastq: UploadFile = File(...),
    wt_fasta: UploadFile = File(...),
    sample_name: str = Form("sample"),
    deletion_start: int = Form(...),
    deletion_end: int = Form(...),
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
    if deletion_start <= 0 or deletion_end <= deletion_start:
        raise HTTPException(400, "Coordenadas inválidas: end debe ser mayor que start (1-based).")

    sample_name = _safe_sample_name(sample_name)
    job_id = datetime.now().strftime("%Y%m%d_%H%M%S_") + uuid.uuid4().hex[:6]
    job_dir = WORK_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    # Guardar uploads
    fastq_name = fastq.filename or "input.fastq.gz"
    if not fastq_name.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
        raise HTTPException(400, "El archivo de lecturas debe ser FASTQ (.fastq/.fastq.gz/.fq/.fq.gz).")
    fastq_path = job_dir / fastq_name
    with open(fastq_path, "wb") as f:
        shutil.copyfileobj(fastq.file, f)

    wt_path = job_dir / "amplicon_WT_input.fa"
    with open(wt_path, "wb") as f:
        shutil.copyfileobj(wt_fasta.file, f)

    params = {
        "deletion_start": deletion_start,
        "deletion_end": deletion_end,
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
        try:
            result = run_full_pipeline(
                job_dir=job_dir,
                fastq_path=fastq_path,
                wt_fasta_path=wt_path,
                sample_name=sample_name,
                params=params,
                log=_log,
                set_progress=_set_progress,
            )
            job["result"] = result
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
            results_dir = job_dir / "results"
            with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
                for p in sorted(results_dir.rglob("*")):
                    if p.is_file():
                        zf.write(p, p.relative_to(results_dir.parent))
        return FileResponse(zip_path, filename=zip_path.name,
                            media_type="application/zip")

    files = job["result"]["files"]
    if kind not in files:
        raise HTTPException(404, f"Tipo de archivo desconocido: {kind}")
    path = Path(files[kind])
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
