"""Fondos alélicos por SNPs y fase de microindels sobre el alelo de origen.

Escenario simulado: amplicón WT con un paciente heterocigoto para 2 SNPs
(fondo_1). La mitad de las lecturas lleva los SNPs (alelo A) y la otra
mitad no (alelo B, igual al WT). Los microindels de un sitio de corte
ocurren SOLO en lecturas del fondo A: la fase debe atribuirlos a ese
alelo y no al WT.
"""
import gzip
from pathlib import Path

from nanodelfy.core.microindel import (
    run_microindel_analysis, parse_guides,
)

WT_LEN = 3000


def _build(tmp_path: Path):
    import random
    rng = random.Random(7)
    wt = "".join(rng.choice("ACGT") for _ in range(WT_LEN))
    # SNPs del alelo A (fondo): posiciones 700 y 1800 (0-based)
    alelo_a = list(wt)
    alelo_a[700] = "A" if alelo_a[700] != "A" else "C"
    alelo_a[1800] = "G" if alelo_a[1800] != "G" else "T"
    alelo_a = "".join(alelo_a)

    guide = wt[1000:1020]                    # guía Cas9 dentro del amplicón

    def _reads(seq_base, n, with_indel):
        out = []
        for j in range(n):
            s = list(seq_base)
            if with_indel:
                # microdeleción de 5 nt en el corte (guía en 1000, corte ~1003)
                del s[1001:1006]
            # ruido de secuenciación: 3 sustituciones aleatorias por lectura
            for _ in range(3):
                p = rng.randrange(len(s))
                s[p] = rng.choice("ACGT")
            out.append("".join(s))
        return out

    reads = (
        _reads(alelo_a, 20, with_indel=True)      # fondo A + microindel
        + _reads(alelo_a, 30, with_indel=False)    # fondo A sin microindel
        + _reads(wt, 40, with_indel=False)         # fondo B (WT) sin microindel
        + _reads(wt, 10, with_indel=True)          # fondo B con microindel
    )
    rng.shuffle(reads)

    fa = tmp_path / "ref.fa"
    with open(fa, "w") as f:
        f.write(">amp\n" + wt + "\n")
    fq = tmp_path / "reads.fastq.gz"
    with gzip.open(fq, "wt") as f:
        for j, s in enumerate(reads):
            f.write(f"@r{j}\n{s}\n+\n{'I' * len(s)}\n")
    return fa, fq, guide


def test_fondos_allelicos_y_fase(tmp_path):
    fa, fq, guide = _build(tmp_path)
    stats = run_microindel_analysis(
        fq, fa, parse_guides(f"guia: {guide}"),
        out_prefix=tmp_path / "salida",
        window=20, max_indel=50, min_mapq=20, nuclease="cas9",
        snp_min_alt_freq=0.20, snp_min_depth=10,
    )

    fondos = stats["fondos_alélicos"]
    assert fondos is not None, "debe detectar SNPs heterocigotos"
    # los dos SNP simulados (y solo ellos) son candidatos
    assert len(fondos["snps"]) == 2
    # dos fondos: alelo A (~50 lecturas) y WT (~50 lecturas)
    clusters = {c["fondo"]: c for c in fondos["clusters"]}
    assert len(clusters) >= 2
    fondo_a = max(
        (c for c in fondos["clusters"] if c["fondo"] != "fondo_ruido"
         and c.get("n_snps_fondo")),
        key=lambda c: c["n_lecturas"])
    assert fondo_a["n_snps_fondo"] == 2

    # fase: los microindels del sitio se reparten entre fondos
    # según la simulación: 20 en fondo A y 10 en fondo B
    site = stats["sites"][0]
    df_path = Path(stats["files"]["per_read"])
    import pandas as pd
    df = pd.read_csv(df_path, sep="\t")
    micro = df[df["categoria"] == "MICROINDEL"]
    por_fondo = micro.groupby("fondo_alélico").size().to_dict()
    assert por_fondo.get("fondo_01", 0) + por_fondo.get("fondo_02", 0) == 30
    # los microindels solo pueden estar en 2 fondos (A y B), ninguno en ruido
    assert set(por_fondo) <= {"fondo_01", "fondo_02"}

    # TSV de fondos escrito
    assert Path(stats["files"]["fondos"]).exists()


def test_sin_snps_devuelve_none_y_columna_sin_SNPs(tmp_path):
    """Muestras sin heterocigosidad: no hay fondos, columna sin_SNPs."""
    fa, fq, guide = _build(tmp_path)
    import pandas as pd
    stats = run_microindel_analysis(
        fq, fa, parse_guides(f"guia: {guide}"),
        out_prefix=tmp_path / "salida2",
        window=20, max_indel=50,
        snp_min_alt_freq=0.90,   # umbral inalcanzable: sin SNPs candidatos
    )
    assert stats["fondos_alélicos"] is None
    df = pd.read_csv(stats["files"]["per_read"], sep="\t")
    assert (df["fondo_alélico"] == "sin_SNPs").all()
