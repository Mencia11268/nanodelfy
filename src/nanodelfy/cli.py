"""Interfaz de línea de comandos de nanodelfy.

Subcomandos:
  nanodelfy run   — ejecuta el pipeline completo sobre una muestra (reproducible)
  nanodelfy web   — arranca la interfaz web local (igual que la app original)

Los parámetros se pasan por opciones CLI o por un JSON de preset
(--preset). Los presets incluidos (nanodelfy/presets/*.json) son ejemplos
reales de los amplicones AESC7 y MmC7 usados en el desarrollo.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Parámetros por defecto: neutros y documentados en el README.
# El usuario SIEMPRE debe indicar las coordenadas de su deleción esperada
# (deletion_start / deletion_end, 1-based inclusivo) en su amplicón.
DEFAULT_PARAMS: dict = {
    "deletion_start": None,        # obligatorio
    "deletion_end": None,          # obligatorio
    "min_quality": 10,
    "min_length": 200,
    "max_length": 5000,
    "breakpoint_tolerance": 50,
    "nhej_search_window": 500,
    "min_junction_insertion": 5,
    "min_sv_deletion": 1500,
    "min_edit_ref_coverage": 0.80,
    "min_mapq": 20,
    "wt_min": None,                # opcional: ventana de longitud de lecturas WT
    "wt_max": None,
    "guides": "",
    "nuclease": "cas9",
    "microindel_window": 20,
    "microindel_max_indel": 50,
    "microindel_large_del_threshold": 100,
    "frt_analysis": None,          # None=auto, True/False explícito
    "frt_guide": None,
}


def _load_preset(path: Path) -> dict:
    with open(path) as f:
        preset = json.load(f)
    if not isinstance(preset, dict):
        raise SystemExit(f"Preset inválido: {path} (debe ser un objeto JSON)")
    return preset


def _cmd_run(args: argparse.Namespace) -> None:
    from .core.pipeline import run_full_pipeline

    params = dict(DEFAULT_PARAMS)
    if args.preset:
        params.update(_load_preset(Path(args.preset)))
    # Las opciones CLI tienen prioridad sobre el preset.
    # Mapeo de nombres argparse → nombres de parámetro del pipeline.
    aliases = {"del_start": "deletion_start", "del_end": "deletion_end"}
    for key, val in vars(args).items():
        if val is None:
            continue
        key = aliases.get(key, key)
        if key in params:
            params[key] = val

    ds, de = params.get("deletion_start"), params.get("deletion_end")
    if ds is None and de is None:
        if not params.get("guides"):
            raise SystemExit(
                "Faltan las coordenadas de la deleción esperada: usa "
                "--del-start/--del-end, un --preset que las incluya o "
                "--guides con dos o más guías para calcularlas "
                "(1-based, inclusivo, sobre el amplicón WT).")
        # ambas vacías con guías: las calcula el pipeline
    elif ds is None or de is None:
        raise SystemExit(
            "Indica ambas coordenadas (--del-start y --del-end) o ninguna "
            "(se calculan desde las guías).")

    # Varios FASTQ en la misma llamada: una muestra por archivo, con las
    # mismas condiciones. El nombre de cada muestra sale de su archivo.
    fastqs = [Path(f) for f in args.fastq]
    if args.name and len(fastqs) > 1:
        raise SystemExit(
            "--name sólo puede usarse con un único FASTQ; con varios, "
            "cada muestra toma el nombre de su archivo.")

    def _nombre(f: Path) -> str:
        if args.name:
            return args.name
        n = f.name
        for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if n.lower().endswith(suf):
                n = n[: -len(suf)]
                break
        return n.split(".")[0] or "muestra"

    try:
        for f in fastqs:
            if len(fastqs) > 1:
                print(f"\n=== Muestra: {f.name} ===")
            res = run_full_pipeline(
                job_dir=Path(args.out) / _nombre(f),
                fastq_path=f,
                wt_fasta_path=Path(args.wt),
                sample_name=_nombre(f),
                params=params,
                log=lambda m: print(m),
            )
    except ValueError as e:
        raise SystemExit(f"Error: {e}")
    if len(fastqs) == 1:
        out_dir = Path(args.out) / _nombre(fastqs[0]) / "results"
        print(f"\nResultados en {out_dir}")
        for k in ("read_classification", "summary", "microindels_summary"):
            if k in res.get("files", {}):
                print(f"  · {res['files'][k]}")
    else:
        print(f"\nResultados en {Path(args.out)} (una subcarpeta por muestra)")


def _cmd_web(args: argparse.Namespace) -> None:
    import threading
    import time
    import webbrowser

    import uvicorn
    from .server import app

    url = f"http://{args.host}:{args.port}"
    print()
    print("=" * 60)
    print(f"  nanodelfy web · abre {url} en tu navegador")
    print("  Pulsa Ctrl+C para detener la aplicación")
    print("=" * 60)
    print()

    if not args.no_browser:
        def _open():
            time.sleep(1.2)
            try:
                webbrowser.open(url)
            except Exception:
                pass
        threading.Thread(target=_open, daemon=True).start()

    uvicorn.run(app, host=args.host, port=args.port, log_level="warning")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="nanodelfy",
        description="nanodelfy — deleciones CRISPR en amplicones de lecturas largas "
                    "(Nanopore/PacBio) con Cas9 o Cas12a.")
    sub = ap.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="ejecutar el pipeline sobre una o varias muestras")
    run.add_argument("--fastq", required=True, nargs="+",
                    help="uno o varios FASTQ (opcionalmente .gz); uno por "
                         "muestra, con las mismas condiciones")
    run.add_argument("--wt", required=True,
                    help="FASTA del amplicón WT (la referencia EDIT se genera sola)")
    run.add_argument("--out", required=True, help="directorio de salida")
    run.add_argument("--name", default=None, help="nombre de la muestra")
    run.add_argument("--preset", default=None,
                    help="JSON con parámetros (ver nanodelfy/presets/)")
    run.add_argument("--del-start", type=int, default=None,
                    help="inicio de la deleción esperada (1-based, inclusivo)")
    run.add_argument("--del-end", type=int, default=None,
                    help="fin de la deleción esperada (1-based, inclusivo)")
    run.add_argument("--min-quality", type=int, default=None)
    run.add_argument("--min-length", type=int, default=None)
    run.add_argument("--max-length", type=int, default=None)
    run.add_argument("--breakpoint-tolerance", type=int, default=None)
    run.add_argument("--min-sv-deletion", type=int, default=None,
                    help="tamaño mínimo de deleción estructural (nt). "
                         "Bájalo (p. ej. 100) para deleciones pequeñas (~120 pb)")
    run.add_argument("--min-mapq", type=int, default=None)
    run.add_argument("--wt-min", type=int, default=None,
                    help="longitud mínima de lecturas 'tipo WT' "
                         "(por defecto: amplicón −4 %%)")
    run.add_argument("--wt-max", type=int, default=None,
                    help="longitud máxima de lecturas 'tipo WT' "
                         "(por defecto: amplicón +4 %%)")
    run.add_argument("--nhej-search-window", type=int, default=None)
    run.add_argument("--guides", default=None,
                    help="guías, una por línea: 'etiqueta: SECUENCIA' "
                         "(ADN o ARN, 20 nt Cas9 / 23-24 nt Cas12a)")
    run.add_argument("--nuclease", choices=["cas9", "cas12a"], default=None,
                    help="nucleasa: cas9 (PAM NGG, corte romo) o cas12a "
                         "(PAM TTTV, corte escalonado 18/23 nt)")
    run.add_argument("--microindel-window", type=int, default=None)
    run.add_argument("--microindel-max-indel", type=int, default=None)
    run.add_argument("--microindel-large-del-threshold", type=int, default=None)
    run.add_argument("--frt-analysis", dest="frt_analysis",
                    action="store_true", default=None,
                    help="forzar el módulo FRT (por defecto: auto)")
    run.add_argument("--no-calibration", dest="calibration",
                    action="store_false", default=None,
                    help="desactivar la capa estadística de calibración "
                         "del ruido (p-valores empíricos con las lecturas WT)")
    run.set_defaults(func=_cmd_run)

    web = sub.add_parser("web", help="interfaz web local")
    web.add_argument("--host", default="127.0.0.1")
    web.add_argument("--port", type=int, default=8765)
    web.add_argument("--no-browser", action="store_true")
    web.set_defaults(func=_cmd_web)

    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
