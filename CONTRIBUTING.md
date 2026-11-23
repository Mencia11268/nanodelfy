# Contribuir a nanodelfy

Gracias por tu interés.

## Desarrollo local

```bash
git clone https://github.com/Mencia11268/nanodelfy.git
cd nanodelfy
pip install -e . pytest
pytest
```

Los tests son autónomos (generan sus propios datos sintéticos) y tardan
~2 segundos. Todo cambio en el clasificador o en los sitios de corte debe
incluir o actualizar tests.

## Qué contribuye

- **Presets**: añade un JSON en `src/nanodelfy/presets/` con los parámetros
  validados para tu amplicón (sin datos personales de pacientes).
- **Nucleasas**: las reglas de corte viven en `_NUCLEASE_SPECS`
  (`src/nanodelfy/core/microindel.py`); añade la tuya con su bibliografía.
- **Bugs**: abre un issue con la versión (`nanodelfy --help` no la muestra;
  usa `python -c "import nanodelfy; print(nanodelfy.__version__)"`), el comando
  exacto y, si es posible, un FASTQ mínimo que reproduzca el fallo.

## Normas

- Código en español o inglés, pero consistente dentro de cada módulo
  (el código existente está comentado en español).
- Mantén la compatibilidad de la CLI y de los presets (1.0.0 es API pública).
- No subas datos de secuenciación reales sin anonimizar.
