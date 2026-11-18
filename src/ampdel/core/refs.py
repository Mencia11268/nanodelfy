"""Genera la referencia editada (post-deleción) a partir del FASTA WT."""
from __future__ import annotations
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def write_fasta_index(fa: Path):
    """Índice FAI en Python puro, sin samtools ni WSL."""
    filas = []
    nombre = None
    longitud = 0
    offset = bases_linea = bytes_linea = 0
    with open(fa, "rb") as f:
        while True:
            linea = f.readline()
            if not linea:
                break
            if linea.startswith(b">"):
                if nombre is not None:
                    filas.append((nombre, longitud, offset, bases_linea, bytes_linea))
                nombre = linea[1:].split()[0].decode("ascii")
                longitud = 0
                offset = f.tell()
                bases_linea = bytes_linea = 0
            else:
                bases = len(linea.rstrip(b"\r\n"))
                if not bases_linea and bases:
                    bases_linea, bytes_linea = bases, len(linea)
                longitud += bases
    if nombre is not None:
        filas.append((nombre, longitud, offset, bases_linea, bytes_linea))
    Path(str(fa) + ".fai").write_text(
        "".join("\t".join(map(str, fila)) + "\n" for fila in filas),
        encoding="ascii")


def build_edited_reference(wt_fa: Path, edit_fa: Path,
                           deletion_start: int, deletion_end: int) -> dict:
    """
    Aplica la deleción [deletion_start, deletion_end] (1-based, inclusiva)
    al amplicón WT y guarda el resultado en edit_fa.
    Devuelve un dict con estadísticas.
    """
    rec = next(SeqIO.parse(str(wt_fa), "fasta"))
    seq = str(rec.seq).upper()
    wt_len = len(seq)

    if not (1 <= deletion_start < deletion_end <= wt_len):
        raise ValueError(
            f"Coordenadas inválidas: start={deletion_start}, end={deletion_end}, "
            f"longitud WT={wt_len}. Debe cumplir 1 ≤ start < end ≤ longitud."
        )

    ds, de = deletion_start, deletion_end
    edit_seq = seq[:ds - 1] + seq[de:]
    expected = wt_len - (de - ds + 1)
    assert len(edit_seq) == expected

    edit_rec = SeqRecord(
        Seq(edit_seq),
        id=f"{rec.id}_EDIT_del{ds}_{de}",
        description=f"WT con deleción de {de - ds + 1} nt entre posiciones {ds}-{de}",
    )
    SeqIO.write([edit_rec], str(edit_fa), "fasta")

    # Indexar ambos FASTA sin programas externos.
    for fa in (wt_fa, edit_fa):
        write_fasta_index(fa)

    return {
        "wt_length": wt_len,
        "deletion_size": de - ds + 1,
        "edit_length": len(edit_seq),
        "wt_id": rec.id,
        "edit_id": edit_rec.id,
    }
