"""Tests de la localización de sitios de corte (Cas9 y Cas12a).

Valores esperados verificados contra el amplicón MmC7 real y la
bibliografía: Cas12a corta con nicks a 18 y 23 nt del extremo
PAM-proximal (Zetsche et al. 2015, Cell 163, 759-771).
"""
import pytest

from nanodelfy.core.microindel import locate_cut_sites, parse_guides, revcomp
from nanodelfy.core.microindel import _normalize_nuclease

# Referencia sintética que reproduce las posiciones del amplicón MmC7:
# PAM TTTC en la posición 46 (1-based) seguida del protospacer Cas12a real.
CAS12A_PROTO = "TATCACTCAGGATGAGGCTGCGA"    # 23 nt, pos 50-72 (1-based)
FRT_GUIDE = "CCGAAGTTCCTATTCTCTAG"           # 20 nt, guía FRT (Cas9-like)
# La guía FRT se coloca en el índice 0-based 1413 para reproducir el corte
# histórico en 1430 (1413 + 20 - 3) del amplicón MmC7 real.
WT = ("A" * 45 + "TTTC" + CAS12A_PROTO + "A" * (1413 - 72)
      + FRT_GUIDE + "A" * 1000)


def test_cas12a_guia_real():
    """Guía Cas12a real: nicks a 18/23 nt, PAM TTTC aguas arriba."""
    sites = locate_cut_sites(WT, parse_guides(f"cas12a_1: {CAS12A_PROTO}"),
                            nuclease="cas12a")
    assert len(sites) == 1
    s = sites[0]
    # protospacer 0-based en 49; nicks 49+18=67 y 49+23=72
    assert s.cut_region == (67, 72)
    assert s.cut_pos == 69                     # centro de la región
    assert s.pam_found and s.pam == "TTTC"
    assert s.nuclease == "cas12a"


def test_cas12a_hebra_negativa():
    """Simétrico en la hebra negativa: PAM aguas abajo en coords +."""
    proto = "GGTTACCACTGATTGACATCTAG"          # 23 nt no palindrómica
    synth = "TTTA" + proto + "TAAA" + "GGCC" * 10
    sites = locate_cut_sites(synth, parse_guides(f"neg: {revcomp(proto)}"),
                            nuclease="cas12a")
    neg = [s for s in sites if s.strand == "-"]
    assert len(neg) == 1
    s = neg[0]
    # p=4, L=23: nicks 4+23-18=9 y 4+23-23=4 → región (4, 9), 5 nt
    assert s.cut_region == (4, 9)
    assert s.pam_found and s.pam == "TTTA"


def test_cas9_regresion():
    """Regla Cas9 sin cambios: corte romo 3 nt aguas arriba del PAM."""
    # guía FRT del MmC7 (20 nt, sin PAM NGG adyacente)
    sites = locate_cut_sites(WT, parse_guides(f"FRT: {FRT_GUIDE}"),
                            nuclease="cas9")
    assert len(sites) == 1
    s = sites[0]
    assert s.cut_region == (s.cut_pos, s.cut_pos)   # romo
    assert not s.pam_found                          # PAM informativo, no bloquea
    assert s.cut_pos == 1430                        # valor verificado pre-cambio


def test_nuclease_por_defecto_es_cas9():
    """Sin parámetro, comportamiento histórico idéntico."""
    sites = locate_cut_sites(WT, parse_guides(f"FRT: {FRT_GUIDE}"))
    assert sites[0].nuclease == "cas9"
    assert sites[0].cut_pos == 1430


def test_alias_y_errores():
    assert _normalize_nuclease("CPF1") == "cas12a"
    assert _normalize_nuclease("AsCas12a") == "cas12a"
    assert _normalize_nuclease("spcas9") == "cas9"
    assert _normalize_nuclease("") == "cas9"
    with pytest.raises(ValueError):
        _normalize_nuclease("cas13")


def test_aviso_protospacer_atipico(capsys):
    """Protospacer de 22 nt en cas12a: aviso, pero distancias desde el PAM."""
    proto22 = "GGTTACCACTGATTGACATCTA"
    synth = "TTTA" + proto22 + "TAAA" + "GGCC" * 10
    sites = locate_cut_sites(synth, parse_guides(f"x: {revcomp(proto22)}"),
                            nuclease="cas12a")
    neg = [s for s in sites if s.strand == "-"][0]
    # p=4, L=22: nicks 4+22-18=8 y 4+22-23=3 → región (3, 8), 5 nt
    assert neg.cut_region == (3, 8)
    assert "AVISO" in capsys.readouterr().out
