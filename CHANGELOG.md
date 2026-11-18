# Changelog

Formato basado en [Keep a Changelog](https://keepachangelog.com/es-ES/1.1.0/).
Versionado [SemVer](https://semver.org/lang/es/).

## [1.0.0] — 2026-11-18

Primera versión pública. ampdel unifica y generaliza la app NanoDel
(versiones AESC7 y MmC7) para su uso con cualquier amplicón.

### Añadido

- Unificación de las dos versiones de NanoDel (AESC7 y MmC7) en un único
  paquete `ampdel` con estructura `src/` y CLI de dos subcomandos:
  `ampdel run` (pipeline reproducible) y `ampdel web` (interfaz web local).
- Generalización: coordenadas de deleción, amplicón y guías son parámetros;
  la referencia EDIT se genera desde el WT automáticamente.
- Parámetro `nuclease` (`cas9`/`cas12a`) en `locate_cut_sites` y en la
  interfaz web: reglas de corte específicas por nucleasa (Cas9: NGG, corte
  romo; Cas12a: TTTV, nicks a 18/23 nt del extremo PAM-proximal).
- Presets JSON (`aesc7.json`, `mmc7.json`) con los parámetros validados de
  los amplicones originales.
- Tests automatizados: sitios de corte (Cas9 y Cas12a, hebras ±) y test de
  humo end-to-end del pipeline con datos sintéticos (30 WT + 30 EDIT,
  deleción de 120 nt).
- Documentación completa: README (español e inglés), CHANGELOG, LICENSE (MIT),
  CITATION.cff, CONTRIBUTING y flujo de CI (GitHub Actions).
- Ventana WT (`wt_min`/`wt_max`) opcional: si no se indica, se deriva del
  tamaño del amplicón (±4 %) en lugar de estar fijada a un amplicón concreto.
- Módulo FRT opcional con activación explícita (`frt_analysis` true/false)
  o automática (guía FRT canónica).

### Corregido

- Parche R3EB_B: reconstrucción de uniones de deleción a partir de colas
  soft-clipped (incluido desde las versiones reparadas de NanoDel).
- Falso positivo en lecturas WT con errores concentrados en un segmento
  (clasificador estructural heredado).

### Heredado de NanoDel (versiones previas)

- Clasificador por lectura con doble alineamiento (WT y EDIT).
- Análisis de microindels alrededor de cortes de guía.
- Módulo FRT específico (MmC7).

## [1.1.0] — 2026-11-18

### Añadido

- Capa estadística de calibración del ruido (no paramétrica, tipo
  CRISPRLungo pero empírica): la distribución nula se estima con las
  propias lecturas WT de la muestra; cada lectura informativa recibe
  p-valores empíricos (`p_as_margin`, `p_max_d`) y un indicador
  `sig_noise`; el resumen incluye CALIB_NULL / CALIB_SIG /
  CALIB_EXP_FP (falsos positivos esperados bajo H0). La clasificación
  determinista NO cambia: la capa es aditiva y se desactiva sola si
  hay menos de 20 lecturas WT (aviso) o con `--no-calibration`.
- Búsqueda directa de la unión (4ª fuente de evidencia) en
  `classify_deletion`: semillas exactas de 8 nt a ambos lados de la
  unión esperada + verificación de identidad de flancos. Resuelve el
  fallback EDIT débil en deleciones pequeñas (~120 pb): la unión
  flanco5+flanco3 es evidencia independiente del tamaño de la deleción.
- Control de especificidad de la búsqueda de unión: si la unión ya
  existe de forma nativa en el WT (flancos dentro de repeticiones
  compuestas, p. ej. sitios FRT pseudopalíndromos), la búsqueda se
  desactiva automáticamente para todo el amplicón (descubierto con el
  control MmC7 real: 8.579 falsos positivos sin este control; 0 con él).

### Validación

- Control MmC7 (9.648 lecturas WT simuladas badread ont_hac): clases
  estructurales idénticas a la versión 1.0.0 (0 falsos positivos);
  capa estadística: 0 lecturas editadas significativas, FP esperados 0.
- amp120 (deleción real de 120 pb, 293 WT + 298 EDIT): sensibilidad
  94,6 % → 95,0 % (una lectura EDIT recuperada por la búsqueda de
  unión), 0 falsos positivos; las 283 lecturas editadas son
  estadísticamente significativas (p < 0,05), FP esperados 1,05.
- Tests: 15/15 (8 de búsqueda de unión y calibración, incluido el caso
  de repeticiones FRT).

### Corregido

- CLI: faltaban las opciones `--wt-min`/`--wt-max` (provocaba un error
  inmediato si se usaban).
