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
