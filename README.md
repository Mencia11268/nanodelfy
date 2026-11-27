# nanodelfy

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24498244.svg)](https://doi.org/10.5281/zenodo.24498244)

**nanodelfy** (de *NanoDel*, su desarrollo original) clasifica lecturas largas
(Nanopore/PacBio) de amplicones CRISPR con deleciones grandes entre dos puntos
de corte, con soporte para **Cas9 y Cas12a**, amplicones de cualquier tamaño y
deleciones desde ~120 pb hasta varios kb.

nanodelfy nace de la unificación de las dos versiones de la app NanoDel
desarrolladas para los amplicones AESC7 y MmC7 (deleciones FRT-FRT), ahora
generalizada para que cualquier laboratorio pueda analizar sus propios
amplicones.

- Dos referencias: solo necesitas el FASTA del amplicón **WT**; nanodelfy genera la
  referencia EDIT (amplicón con la deleción esperada) automáticamente.
- Clasificación por lectura: WT, deleción exacta, deleción extendida, deleción con
  inserción, inversión, lecturas truncadas (no informativas)...
- Microindels alrededor de los cortes de las guías, con reglas específicas de
  Cas9 (PAM NGG, corte romo) y Cas12a (PAM TTTV, corte escalonado).
- Dos interfaces: **web local** (arrastra y suelta tu FASTQ) y **CLI**
  reproducible (`nanodelfy run`), con presets JSON.

## Versión web alojada (sin instalación)

Puedes usar la interfaz web de nanodelfy sin instalar nada, desplegándola
gratis en tu propia cuenta de Render (no requiere tarjeta; el servicio se
duerme tras 15 minutos de inactividad y despierta solo al recibir una
visita):

[![Deploy to Render](https://render.com/images/deploy-to-render-button.svg)](https://render.com/deploy?repo=https://github.com/Mencia11268/nanodelfy)

Pulsa el botón, crea una cuenta gratuita (o inicia sesión), confirma el
plan gratuito y espera unos minutos a que se construya. Obtendrás una URL
pública tipo `https://nanodelfy-xxxx.onrender.com` lista para compartir con
otros usuarios. Sube las muestras en formato `.fastq.gz` (comprimido) para
minimizar el tiempo de subida.

## Instalación

Requisitos: Python 3.10 o superior.

```bash
pip install nanodelfy          # desde PyPI (cuando se publique)
```

Desde el código fuente (este repositorio):

```bash
git clone https://github.com/Mencia11268/nanodelfy.git
cd nanodelfy
pip install .
```

Dependencias: `mappy` (alineador minimap2 para Python), `biopython`, `pandas`,
`openpyxl`, `fastapi`, `uvicorn`. La escritura de BAM es en Python puro (no
necesita samtools ni pysam).

En **Windows**, `mappy` no tiene binarios precompilados: usa el subsistema WSL
(o conda) para instalar nanodelfy.

## Inicio rápido

### Interfaz web

```bash
nanodelfy web            # abre http://127.0.0.1:8765
```

Sube tu FASTQ y el FASTA del amplicón WT, indica las coordenadas de la deleción
esperada y las guías, y pulsa Ejecutar.

### Línea de comandos

```bash
nanodelfy run \
  --fastq muestra.fastq.gz \
  --wt amplicon_WT.fa \
  --out salida \
  --del-start 1438 --del-end 1557 \
  --min-sv-deletion 100 \
  --guides "guia1: CCGAAGTTCCTATTCTCTAG" \
  --nuclease cas9
```

O con un preset:

```bash
nanodelfy run --fastq muestra.fastq --wt amplicon_WT.fa --out salida \
  --preset $(python -c "import nanodelfy, pathlib; print(pathlib.Path(nanodelfy.__file__).parent / 'presets' / 'mmc7.json')")
```

## Parámetros

Los mismos parámetros funcionan en la web (campos del formulario) y en el CLI.

| Parámetro | Por defecto | Descripción |
|---|---|---|
| `deletion_start`, `deletion_end` | obligatorio | Coordenadas 1-based inclusivas de la deleción esperada sobre el amplicón WT. **Es el único dato obligatorio.** |
| `min_quality` | 10 | Calidad media mínima por lectura (Phred). Bájalo a 7 para datos ruidosos. |
| `min_length`, `max_length` | 200 / 5000 | Longitud mínima/máxima de lectura tras el filtrado. Ajusta al rango esperado de tu amplicón (p. ej. 800-5000 para amplicones de ~3,9 kb). |
| `breakpoint_tolerance` | 50 | Tolerancia (nt) para considerar una unión de deleción "exacta". |
| `min_sv_deletion` | 1500 | Tamaño mínimo (nt) de la deleción estructural detectable. **Bájalo por debajo del tamaño de tu deleción**: 100 para deleciones de ~120 pb. |
| `min_mapq` | 20 | MAPQ mínimo del alineamiento primario. |
| `min_edit_ref_coverage` | 0.80 | Cobertura mínima de la referencia EDIT para aceptar una lectura como editada. |
| `wt_min`, `wt_max` | auto (±4 %) | Ventana de longitud de lecturas "tipo WT". Si se omite, se deriva del tamaño del amplicón. |
| `nhej_search_window` | 500 | Ventana de búsqueda de uniones NHEJ alrededor de cada breakpoint. |
| `guides` | vacío | Guías, una por línea en formato `etiqueta: SECUENCIA` (ADN o ARN, U/T). Con dos o más guías se calcula además el tamaño teórico de la deleción (distancia entre los cortes, según la nucleasa), que se muestra en el gráfico de uniones junto a la diferencia (Δ) de cada variante observada. |
| `nuclease` | cas9 | `cas9` (protospacer 20 nt, PAM NGG, corte romo 3 nt aguas arriba del PAM) o `cas12a` (protospacer 23-24 nt, PAM TTTV, nicks a 18 y 23 nt del extremo PAM-proximal). |
| `microindel_window` | 20 | Ventana (nt) alrededor del corte para el análisis de microindels. |
| `microindel_max_indel` | 50 | Indel máximo (nt) considerado microindel. |
| `microindel_large_del_threshold` | 100 | Umbral (nt) entre microdeleción y deleción grande en las tablas de microindels. |
| `frt_analysis` | auto | Módulo FRT opcional (análisis conservador de la región FRT, herencia MmC7). Auto = se activa si una guía es la FRT canónica. |

### Valores recomendados según el tamaño de la deleción

| Tamaño de deleción | `min_sv_deletion` | `breakpoint_tolerance` | Notas |
|---|---|---|---|
| Grande (>1 kb, p. ej. FRT-FRT) | 1500 | 50 | Valores por defecto; presets AESC7/MmC7. |
| Media (0,3-1 kb) | ~200 | 30 | Sube `min_mapq` si hay falsos positivos. |
| Pequeña (~120 pb) | 100 | 20 | Amplía `nhej_search_window` a 200. |

## Salidas

En `salida/results/`:

- `<muestra>.read_classification.tsv` — clase por lectura: `read_id`, `length`,
  `class`, alineamientos WT/EDIT, scores, regla aplicada, offsets respecto a los
  breakpoints, inserciones; desde la v1.1.0, columnas de calibración
  (`max_d_wt`, `as_margin`, `p_as_margin`, `p_max_d`, `sig_noise`).
- `<muestra>.summary.tsv` — recuentos y porcentajes por clase (incluye el
  denominador de lecturas informativas).
- `<muestra>.resultados.xlsx` — todo lo anterior en Excel.
- `<muestra>.WT.bam`, `<muestra>.EDIT.bam` (+ índices) — alineamientos por
  referencia, listos para IGV.
- `<muestra>.microindels_*.tsv` — microindels por lectura, resumen por guía y
  alelos alrededor de los cortes.
- `<muestra>.BAMs_para_IGV.zip` — BAMs + FASTAs para abrir en IGV.

### Clases de lectura

| Clase | Significado |
|---|---|
| `WT` / `WT_PARTIAL` | Lectura sin edición (la parcial atraviesa un breakpoint sin deleción; informativa). |
| `DEL_EXACT` | Deleción con unión dentro de la tolerancia de los breakpoints esperados. |
| `DEL_EXTENDED` | Deleción con unión extendida más allá de los breakpoints (dentro de `nhej_search_window`). |
| `DEL_WITH_INSERTION` | Deleción con inserción ≥ `min_junction_insertion` nt en la unión (NHEJ). |
| `DEL_TRUNCATED` | Deleción en lectura truncada. |
| `INV` | Inversión del fragmento entre los breakpoints. |
| `WT_TRUNCATED`, `EDIT_TRUNCATED`, `AMBIGUOUS_TRUNCATED` | Lecturas truncadas no informativas (fuera del denominador). |
| `OTHER` | Indels imprecisos / no informativas. |

## Nucleasas: Cas9 y Cas12a

`locate_cut_sites(..., nuclease=...)` calcula los sitios de corte según la nucleasa:

- **Cas9**: protospacer de 20 nt; corte romo 3 nt aguas arriba del PAM (NGG).
  El parámetro PAM es informativo (no bloquea el análisis si falta).
- **Cas12a** (Cpf1): protospacer de 23-24 nt; PAM TTTV aguas arriba (5');
  corte escalonado con nicks a 18 y 23 nt del extremo PAM-proximal
  (Zetsche et al. 2015). Se aceptan los alias `cpf1`, `ascas12a`, `lbcas12a`.

La búsqueda de la guía es exacta en ambas hebras del amplicón WT; las distancias
al PAM se calculan desde el extremo PAM-proximal del protospacer, por lo que
funciona con guías de 23 o 24 nt.

## Presets incluidos

En `src/nanodelfy/presets/`:

- `aesc7.json` — amplicón AESC7 (3117 nt), deleción FRT-FRT 1157-2981.
- `mmc7.json` — amplicón MmC7 (3887 nt), deleción FRT-FRT 1438-3295 (activa el
  módulo FRT).

Crea los tuyos copiando la estructura (JSON plano con los parámetros de la tabla
anterior) y pásalos con `--preset mi_preset.json`.

## Capa estadística (calibración del ruido)

Desde la v1.1.0 nanodelfy separa señal de ruido **sin alterar la clasificación
determinista**. El método es una calibración empírica:

1. Hipótesis nula H0: "la lectura es WT (solo ruido de secuenciación)".
2. La distribución nula se estima con las **propias lecturas WT** de la muestra
   (clases WT y WT_PARTIAL): comparten el mismo ruido de secuenciación.
   Si hay menos de 20 lecturas WT, la capa se omite con un aviso (o
   desactívala con `--no-calibration`).
3. Estadísticos por lectura: `max_d_wt` (la mayor deleción del CIGAR contra
   WT) y `as_margin` (AS de EDIT − AS de WT).
4. Cada lectura informativa recibe p-valores empíricos un lado
   (`p_as_margin`, `p_max_d`) y el indicador `sig_noise` (p < 0,05).

En el resumen, las filas `CALIB_NULL` (lecturas WT usadas como nula),
`CALIB_SIG` (lecturas editadas con p < 0,05) y `CALIB_EXP_FP` (falsos
positivos esperados bajo H0: la suma de p-valores de las lecturas
editadas — cuántas de ellas serían compatibles con puro ruido WT).

Las lecturas truncadas en regiones idénticas WT/EDIT tienen margen ~0 por
construcción: sus columnas quedan vacías y no cuentan como significativas.

Validación: en una muestra WT pura (9.648 lecturas simuladas), 0 lecturas
significativas y FP esperados 0,0; en la mezcla amp120 (deleción real de
120 pb), las 283 lecturas editadas son significativas con 1,05 FP esperados.

## Búsqueda directa de la unión (deleciones pequeñas)

En deleciones pequeñas (<150 pb) las referencias WT y EDIT son casi idénticas
y el margen de score apenas supera el ruido. Desde la v1.1.0 hay una 4ª fuente
de evidencia: la **unión** (flanco 5' + flanco 3' pegados) se busca
directamente en la secuencia de la lectura con semillas exactas de 8 nt y
verificación de identidad de los flancos. La unión solo existe en lecturas
editadas, sea cual sea el tamaño de la deleción.

Control de especificidad: si la unión ya existe de forma nativa en el
amplicón WT (flancos dentro de repeticiones compuestas, como los sitios FRT
pseudopalíndromos), la señal no es evidencia y la búsqueda se desactiva
automáticamente para todo el amplicón.

## Limitaciones

- La calibración estadística es empírica (no paramétrica): necesita ≥20
  lecturas WT en la muestra; si la edición es casi total (cero WT),
  usa una muestra control WT o desactiva la capa.
- La referencia EDIT se genera como deleción simple: para alelos con inserciones
  grandes conocidas, esas lecturas caerán en `DEL_WITH_INSERTION`/`OTHER`.
- Un solo amplicón por ejecución (por diseño: una pareja de guías y una
  deleción esperada).

## Desarrollo

```bash
pip install -e . pytest
pytest                    # 7 tests: sitios de corte (Cas9/Cas12a) + pipeline end-to-end
```

Los tests usan datos sintéticos autogenerados (no requieren descargas).

## Citación

Si usas nanodelfy, cita la versión archivada en Zenodo:

> Mencia, A., Murillas, R., & Garcia Diez, M. (2026). nanodelfy 1.1.3.
> Zenodo. https://doi.org/10.5281/zenodo.24534251

Para referirse al conjunto de versiones, usa el DOI de concepto:
https://doi.org/10.5281/zenodo.24498243. El archivo `CITATION.cff` contiene los
metadatos de citación.

## Licencia

MIT — ver [LICENSE](LICENSE).

## Historia

nanodelfy unifica y generaliza la app NanoDel (versiones AESC7 y MmC7).
El desarrollo previo está documentado en [CHANGELOG.md](CHANGELOG.md).
