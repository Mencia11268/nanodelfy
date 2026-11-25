"""Regresión del JavaScript de la página de resultados (job.html).

La fila CALIB_EXP_FP del resumen de clasificación llega con
frequency_pct vacío (''); un ''.toFixed(2) rompía el renderizado de la
tabla Y del gráfico de barras de clases. Este test ejecuta el script de
la página en un DOM simulado (node) con un JSON realista y comprueba que
no hay errores y que ambos gráficos se pintan.

Requiere node en el sistema; se omite si no está disponible.
"""
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
JOB_HTML = REPO / "src" / "nanodelfy" / "templates" / "job.html"

HARNESS = r"""
const fs = require('fs');
const jobJson = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));

function makeEl(id) {
  return {
    id, _html: '', style: {}, href: '', textContent: '',
    set innerHTML(v) { this._html = v; }, get innerHTML() { return this._html; },
    appendChild(c) { this._html += (c && c._html) || ''; },
  };
}
const elements = {};
global.document = {
  getElementById: (id) => (elements[id] = elements[id] || makeEl(id)),
  querySelector: () => makeEl('tbody'),
  createElement: () => makeEl('tr'),
  addEventListener: () => {},
};
global.location = { origin: 'http://test', pathname: '/job/x' };
global.fetch = async () => ({ ok: true, json: async () => jobJson });
const errores = [];
console.error = (...a) => errores.push(String(a[1] && a[1].message || a));

eval(fs.readFileSync(process.argv[3], 'utf8'));

setTimeout(() => {
  const g = id => (elements[id] || {})._html || '';
  const cls = g('chart-classes');
  const mic = g('chart-micro');
  const dj = g('deljunc-wrap');
  const iseq = g('indelseq-wrap');
  console.log('RESULT ' + JSON.stringify({
    errores,
    clases_svg: cls.includes('<svg'),
    clases_barras: (cls.match(/<rect/g) || []).length,
    micro_svg: mic.includes('<svg'),
    deljunc_bloques: (dj.match(/aln-block/g) || []).length,
    deljunc_corte: (dj.match(/aln-cut/g) || []).length,
    deljunc_del: (dj.match(/aln-d/g) || []).length,
    indelseq_bloques: (iseq.match(/aln-block/g) || []).length,
    indelseq_del: (iseq.match(/aln-d/g) || []).length,
    indelseq_ins: (iseq.match(/aln-i/g) || []).length,
    indelseq_corte: (iseq.match(/aln-cut/g) || []).length,
  }));
  process.exit(0);
}, 300);
"""

JOB = {
    "status": "done",
    "error": None,
    "params": {},
    "result": {
        "classification": {
            "summary": [
                {"class": "WT", "description": "salvaje", "reads": 50,
                 "frequency_pct": 50.0},
                {"class": "CALIB_EXP_FP",
                 "description": "Falsos positivos esperados bajo H0",
                 "reads": 0.99, "frequency_pct": ""},
                {"class": "TOTAL", "description": "total", "reads": 100,
                 "frequency_pct": 100.0},
            ],
        },
        "microindels": {
            "summary": [
                {"sitio": "guia_1@2017", "con_microindel": 3,
                 "pct_microindel_sobre_evaluables": "4.00"},
            ],
            "alleles_display": [
                {
                    "sitio": "guia_1@2017",
                    "corte_1based": "entre 2017 y 2018",
                    "ref": "ACGTACGTACGTACGTACGT",
                    "corte_rel": 10,
                    "alelos": [
                        {"alelo": "alelo_001", "n_lecturas": 5, "pct": "6.25",
                         "eventos": "DEL3nt@ref2014",
                         "segmentos": [
                             {"t": "m", "seq": "ACGTACGTAC"},
                             {"t": "d", "seq": "GTA"},
                             {"t": "i", "seq": "TT"},
                             {"t": "m", "seq": "CGTACGTACGT"},
                         ]},
                    ],
                },
            ],
        },
        "deletion_view": {
            "ref_len": 3000,
            "esperada": {"start": 999, "end": 2000},
            "cortes": [{"sitio": "guia_1@1005", "pos": 1005}],
            "deleciones": [
                {"start": 1000, "end": 1999, "length": 999, "n": 40,
                 "pct": 40.0, "ins_total": 0, "desde_guia": "guia_1@1005",
                 "hasta_guia": "", "etiqueta": "guia_1@1005",
                 "flanco5": "AAAACCCCGGGGTTTTAAAACCCCGGGGT",
                 "flanco3": "TTTTGGGGCCCCAAAATTTTGGGGCCCCA",
                 "deletado": ""},
            ],
        },
        "files": {},
    },
}


@pytest.mark.skipif(shutil.which("node") is None, reason="node no disponible")
def test_job_page_js_renderiza_graficos(tmp_path):
    html = JOB_HTML.read_text()
    scripts = re.findall(r"<script>(.*?)</script>", html, re.S)
    assert scripts, "job.html debe contener un script"
    script = tmp_path / "job.js"
    script.write_text(scripts[0])
    harness = tmp_path / "harness.js"
    harness.write_text(HARNESS)
    job_json = tmp_path / "job.json"
    job_json.write_text(json.dumps(JOB))
    out = subprocess.run(
        ["node", str(harness), str(job_json), str(script)],
        capture_output=True, text=True, timeout=60, check=True,
    )
    result = json.loads(
        [l for l in out.stdout.splitlines() if l.startswith("RESULT ")][0][7:]
    )
    assert result["errores"] == [], f"errores JS: {result['errores']}"
    assert result["clases_svg"], "el gráfico de clases debe pintarse"
    assert result["clases_barras"] >= 2, "deben pintarse barras (WT + CALIB)"
    assert result["micro_svg"], "el gráfico de microindels debe pintarse"
    # secuencia de la unión de deleciones recurrentes
    assert result["deljunc_bloques"] == 1, "debe haber 1 bloque de unión"
    assert result["deljunc_corte"] == 2, \
        "las dos bases de la unión deben resaltarse"
    assert result["deljunc_del"] >= 1, "el segmento eliminado debe ir en rojo"
    # alineamiento de indels recurrentes por guía
    assert result["indelseq_bloques"] == 1, "debe haber 1 bloque de guía"
    assert result["indelseq_del"] >= 1, "la deleción debe ir en rojo"
    assert result["indelseq_ins"] >= 1, "la inserción debe ir en azul"
    assert result["indelseq_corte"] == 1, "el corte debe resaltarse"
