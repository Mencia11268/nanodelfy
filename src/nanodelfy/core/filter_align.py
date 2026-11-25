"""Filtrado y alineamiento Windows-friendly (sin pysam ni samtools).

Usa:
  - mappy   → bindings de minimap2 (instalable con pip en Windows)
  - escritor BAM propio en Python puro (bgzf + estructura BAM) para escritura

Devuelve dos cosas:
  - lista de objetos AlnRecord (con API mínima compatible con pysam.AlignedSegment)
  - ficheros .bam ordenados + índice .bai en disco para IGV
"""
from __future__ import annotations
import gzip, math, struct, zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple, Optional

import mappy as mp


# ----------------------------- AlnRecord -----------------------------
# Imita la subconjunto de la API de pysam.AlignedSegment que necesita
# classify.py. Permite eliminar la dependencia de pysam.

@dataclass
class AlnRecord:
    query_name: str
    flag: int
    reference_name: str
    reference_start: int
    mapping_quality: int
    cigartuples: list           # list[(op, length)]
    query_length: int
    _tags: dict = field(default_factory=dict)
    is_unmapped_flag: bool = False
    # Secuencia de la lectura EN ORIENTACIÓN DE REFERENCIA (igual que el
    # campo SEQ del BAM): para lecturas en hebra reversa es el revcomp de la
    # secuencia original. Los soft-clips del CIGAR se indexan sobre ella.
    # Se usa para reconstruir uniones de deleción a partir de colas
    # soft-clipped (parche R3EB_B).
    query_sequence: Optional[str] = None

    # ---- flags ----
    @property
    def is_reverse(self) -> bool: return bool(self.flag & 16)
    @property
    def is_secondary(self) -> bool: return bool(self.flag & 256)
    @property
    def is_supplementary(self) -> bool: return bool(self.flag & 2048)
    @property
    def is_unmapped(self) -> bool: return bool(self.flag & 4) or self.is_unmapped_flag

    # ---- coordenadas ----
    @property
    def reference_end(self) -> int:
        end = self.reference_start
        for op, L in self.cigartuples or []:
            if op in (0, 2, 3, 7, 8):
                end += L
        return end

    @property
    def reference_length(self) -> int:
        return self.reference_end - self.reference_start

    @property
    def cigarstring(self) -> str:
        codes = "MIDNSHP=X"
        return "".join(f"{L}{codes[op]}" for op, L in (self.cigartuples or []))

    def infer_read_length(self) -> int:
        if not self.cigartuples:
            return self.query_length or 0
        return sum(L for op, L in self.cigartuples if op in (0, 1, 4, 7, 8))

    # ---- coordenadas sobre la lectura (query) ----
    # Imitan pysam.AlignedSegment.query_alignment_start / _end:
    # posición en la lectura donde EMPIEZA y TERMINA la región alineada
    # (excluye los soft/hard clips inicial y final).
    @property
    def query_alignment_start(self) -> int:
        start = 0
        for op, L in self.cigartuples or []:
            if op in (4, 5):   # S o H al inicio
                start += L
            else:
                break
        return start

    @property
    def query_alignment_end(self) -> int:
        qlen = self.infer_read_length()
        trailing = 0
        for op, L in reversed(self.cigartuples or []):
            if op in (4, 5):
                trailing += L
            else:
                break
        return qlen - trailing

    def get_tag(self, tag: str):
        if tag not in self._tags:
            raise KeyError(tag)
        return self._tags[tag]

    def has_tag(self, tag: str) -> bool:
        return tag in self._tags


# ----------------------------- FILTRADO -----------------------------

def filter_reads(in_fastq: Path, out_fastq_gz: Path,
                 min_quality: int, min_length: int, max_length: int,
                 threads: int = 4, log = print) -> dict:
    log(f"  · filtrando Q≥{min_quality}, longitud {min_length}-{max_length} nt")
    n_in, n_out = _filter_python(in_fastq, out_fastq_gz,
                                 min_quality, min_length, max_length)
    log(f"  · {n_in} → {n_out} lecturas tras filtrado")
    return {"reads_in": n_in, "reads_out": n_out}


def _mean_qscore(qual_str: str) -> float:
    if not qual_str:
        return 0.0
    probs = [10 ** (-(ord(c) - 33) / 10.0) for c in qual_str]
    mean_p = sum(probs) / len(probs)
    if mean_p <= 0:
        return 60.0
    return -10.0 * math.log10(mean_p)


def _filter_python(in_fastq, out_fastq_gz, minq, minl, maxl):
    opener = gzip.open if str(in_fastq).endswith(".gz") else open
    n_in = n_out = 0
    with opener(in_fastq, "rt") as fh, gzip.open(out_fastq_gz, "wt") as out:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                break
            n_in += 1
            seq_strip = seq.rstrip("\n")
            qual_strip = qual.rstrip("\n")
            L = len(seq_strip)
            if L < minl or L > maxl:
                continue
            if _mean_qscore(qual_strip) < minq:
                continue
            out.write(header); out.write(seq); out.write(plus); out.write(qual)
            n_out += 1
    return n_in, n_out


# ----------------------------- ALINEAMIENTO -----------------------------

def align_minimap2(fastq_gz: Path, ref_fa: Path, out_bam: Path,
                   threads: int = 4, log = print) -> List[AlnRecord]:
    """Alinea con mappy y escribe un BAM ordenado con índice .bai (Python puro).

    Devuelve la lista de AlnRecord (con API tipo pysam.AlignedSegment) para que
    classify.py pueda usarlos directamente sin tener que releer el BAM.
    """
    log(f"  · alineando contra {ref_fa.name}")

    # bw por defecto de minimap2: fijar bw/bw_long a 2000 rompía el encadenado
    # a través de deleciones grandes (la lectura quedaba parcial: 1438M + clip)
    aligner = mp.Aligner(
        str(ref_fa), preset="map-ont",
        n_threads=threads, best_n=5,
    )
    if not aligner:
        raise RuntimeError(f"No se pudo construir el índice de {ref_fa}")

    ref_lengths = [(name, len(aligner.seq(name))) for name in aligner.seq_names]
    ref_name_to_tid = {n: i for i, (n, _) in enumerate(ref_lengths)}

    # Recolectar alineamientos en memoria (amplicones pequeños)
    alignments: List[dict] = []
    records: List[AlnRecord] = []
    n_reads = 0
    for name, seq, qual in mp.fastx_read(str(fastq_gz), read_comment=False):
        n_reads += 1
        hits = list(aligner.map(seq, cs=False, MD=True))
        if not hits:
            alignments.append({
                "qname": name, "flag": 4, "tid": -1, "pos": -1,
                "mapq": 0, "cigar": [], "seq": seq, "qual": qual,
                "tags": [], "ref_len": 0,
            })
            records.append(AlnRecord(
                query_name=name, flag=4, reference_name="",
                reference_start=-1, mapping_quality=0,
                cigartuples=[], query_length=len(seq),
                _tags={}, is_unmapped_flag=True,
                query_sequence=seq,
            ))
            continue
        primary_seen = False
        for h in hits:
            is_primary = h.is_primary and not primary_seen
            if is_primary:
                primary_seen = True
            flag = 0
            if h.strand == -1:
                flag |= 16
            if not is_primary:
                flag |= 2048  # supplementary
            # CIGAR: mappy entrega [length, op]; añadir soft-clips
            qlen = len(seq)
            if h.strand == 1:
                left_clip = h.q_st
                right_clip = qlen - h.q_en
                bam_seq = seq
                bam_qual = qual
            else:
                left_clip = qlen - h.q_en
                right_clip = h.q_st
                bam_seq = mp.revcomp(seq)
                bam_qual = qual[::-1] if qual else qual
            cigar = []
            if left_clip > 0:
                cigar.append((4, left_clip))
            cigar.extend((op, length) for length, op in h.cigar)
            if right_clip > 0:
                cigar.append((4, right_clip))

            # Calcular ref_len a partir del CIGAR
            ref_len = sum(L for op, L in cigar if op in (0, 2, 3, 7, 8))

            tags = [("NM", "i", h.NM), ("AS", "i", h.mlen)]
            if getattr(h, "MD", None):
                tags.append(("MD", "Z", h.MD))

            alignments.append({
                "qname": name, "flag": flag,
                "tid": ref_name_to_tid[h.ctg],
                "pos": h.r_st,
                "mapq": h.mapq,
                "cigar": cigar,
                "seq": bam_seq, "qual": bam_qual,
                "tags": tags,
                "ref_len": ref_len,
            })
            records.append(AlnRecord(
                query_name=name, flag=flag,
                reference_name=h.ctg,
                reference_start=h.r_st,
                mapping_quality=h.mapq,
                cigartuples=cigar,
                query_length=len(seq),
                _tags={t: v for t, _ttype, v in tags},
                query_sequence=bam_seq,
            ))

    # Ordenar por (tid, pos); las no mapeadas al final
    def sort_key(a):
        if a["tid"] < 0:
            return (1 << 30, 0)
        return (a["tid"], a["pos"])
    alignments.sort(key=sort_key)

    # Escribir BAM
    _write_bam(out_bam, ref_lengths, alignments)
    # Escribir índice .bai
    _write_bai(out_bam, ref_lengths, alignments)
    log(f"  · {n_reads} lecturas → {len(alignments)} alineamientos ({out_bam.name})")
    return records


# ----------------------------- BAM WRITER (Python puro) -----------------------------
#
# Especificación SAM/BAM: https://samtools.github.io/hts-specs/SAMv1.pdf
# Formato BGZF: bloques gzip de máx 65536 bytes, con extra field "BC"
# que indica el tamaño del bloque comprimido.

BGZF_MAX_BLOCK = 65280  # tamaño de datos por bloque (margen de seguridad)
BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)
BAM_MAGIC = b"BAM\x01"

# Codificación BAM de bases: =ACMGRSVTWYHKDBN → 0..15
_BAM_SEQ_TABLE = {c: i for i, c in enumerate("=ACMGRSVTWYHKDBN")}


def _bgzf_compress_block(data: bytes) -> bytes:
    """Comprime un bloque y devuelve el bloque BGZF completo."""
    compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
    cdata = compressor.compress(data) + compressor.flush()
    crc = zlib.crc32(data) & 0xFFFFFFFF
    isize = len(data) & 0xFFFFFFFF
    # BSIZE = tamaño total del bloque - 1
    # Total = 12 (header fijo) + 6 (subfield BC) + len(cdata) + 8 (CRC+ISIZE)
    bsize = 12 + 6 + len(cdata) + 8 - 1
    header = (
        b"\x1f\x8b\x08\x04"                # magic + FLG (FEXTRA)
        + b"\x00\x00\x00\x00"              # MTIME
        + b"\x00"                          # XFL
        + b"\xff"                          # OS
        + struct.pack("<H", 6)             # XLEN = 6
        + b"BC"                            # SI1 SI2
        + struct.pack("<H", 2)             # SLEN = 2
        + struct.pack("<H", bsize)         # BSIZE
    )
    footer = struct.pack("<II", crc, isize)
    return header + cdata + footer


class _BgzfWriter:
    """Escritor BGZF con tracking de virtual offsets para .bai."""

    def __init__(self, path: Path):
        self.fh = open(path, "wb")
        self.buf = bytearray()
        self.coffset = 0  # offset comprimido del bloque actual

    def vo(self) -> int:
        """Virtual offset actual = coffset << 16 | uoffset_in_block."""
        return (self.coffset << 16) | len(self.buf)

    def write(self, data: bytes):
        i = 0
        while i < len(data):
            space = BGZF_MAX_BLOCK - len(self.buf)
            chunk = data[i:i + space]
            self.buf.extend(chunk)
            i += len(chunk)
            if len(self.buf) >= BGZF_MAX_BLOCK:
                self._flush_block()

    def _flush_block(self):
        if not self.buf:
            return
        block = _bgzf_compress_block(bytes(self.buf))
        self.fh.write(block)
        self.coffset += len(block)
        self.buf.clear()

    def close(self):
        self._flush_block()
        self.fh.write(BGZF_EOF)
        self.fh.close()


def _encode_seq(seq: str) -> bytes:
    seq = seq.upper()
    out = bytearray()
    for i in range(0, len(seq), 2):
        a = _BAM_SEQ_TABLE.get(seq[i], 15)
        b = _BAM_SEQ_TABLE.get(seq[i+1], 0) if i + 1 < len(seq) else 0
        out.append((a << 4) | b)
    return bytes(out)


def _encode_qual(qual: str, seq_len: int) -> bytes:
    if not qual:
        return bytes([0xFF] * seq_len)
    return bytes(ord(c) - 33 for c in qual)


def _encode_cigar(cigar: list) -> bytes:
    """CIGAR en BAM: array uint32, ((length << 4) | op)."""
    return b"".join(struct.pack("<I", (L << 4) | op) for op, L in cigar)


def _encode_tags(tags: list) -> bytes:
    out = bytearray()
    for tag, ttype, val in tags:
        out.extend(tag.encode("ascii"))
        if ttype == "i":
            v = int(val)
            if -128 <= v <= 127:
                out.extend(b"c"); out.extend(struct.pack("<b", v))
            elif -32768 <= v <= 32767:
                out.extend(b"s"); out.extend(struct.pack("<h", v))
            else:
                out.extend(b"i"); out.extend(struct.pack("<i", v))
        elif ttype == "f":
            out.extend(b"f"); out.extend(struct.pack("<f", float(val)))
        elif ttype == "Z":
            out.extend(b"Z")
            out.extend(str(val).encode("ascii", errors="replace"))
            out.append(0)
        else:
            raise ValueError(f"tag type no soportado: {ttype}")
    return bytes(out)


def _write_bam(path: Path, ref_lengths, alignments):
    w = _BgzfWriter(path)
    # Header
    hdr_text = "@HD\tVN:1.6\tSO:coordinate\n"
    for name, L in ref_lengths:
        hdr_text += f"@SQ\tSN:{name}\tLN:{L}\n"
    hdr_text += "@PG\tID:nanodel\tPN:nanodel\tCL:mappy minimap2\n"
    hdr_bytes = hdr_text.encode("ascii")

    w.write(BAM_MAGIC)
    w.write(struct.pack("<i", len(hdr_bytes)))
    w.write(hdr_bytes)
    w.write(struct.pack("<i", len(ref_lengths)))
    for name, L in ref_lengths:
        name_b = name.encode("ascii") + b"\x00"
        w.write(struct.pack("<i", len(name_b)))
        w.write(name_b)
        w.write(struct.pack("<i", L))

    # Alineamientos
    for a in alignments:
        qname_b = a["qname"].encode("ascii") + b"\x00"
        cigar_b = _encode_cigar(a["cigar"])
        seq_b = _encode_seq(a["seq"])
        qual_b = _encode_qual(a["qual"] or "", len(a["seq"]))
        tags_b = _encode_tags(a["tags"])

        bin_val = _reg2bin(a["pos"], a["pos"] + max(1, a["ref_len"])) if a["tid"] >= 0 else 4680
        bin_mq_nl = (bin_val << 16) | (a["mapq"] << 8) | len(qname_b)
        flag_nc = (a["flag"] << 16) | len(a["cigar"])

        body = (
            struct.pack("<i", a["tid"])
            + struct.pack("<i", a["pos"])
            + struct.pack("<I", bin_mq_nl)
            + struct.pack("<I", flag_nc)
            + struct.pack("<i", len(a["seq"]))
            + struct.pack("<i", -1)    # next_ref_id
            + struct.pack("<i", -1)    # next_pos
            + struct.pack("<i", 0)     # tlen
            + qname_b
            + cigar_b
            + seq_b
            + qual_b
            + tags_b
        )
        # Record virtual offset (inicio del registro) — útil para .bai
        a["_vo_start"] = w.vo()
        w.write(struct.pack("<i", len(body)))
        w.write(body)
        a["_vo_end"] = w.vo()

    w.close()


def _reg2bin(beg: int, end: int) -> int:
    """Función oficial reg2bin del SAM/BAM spec."""
    end -= 1
    if (beg >> 14) == (end >> 14): return ((1 << 15) - 1) // 7 + (beg >> 14)
    if (beg >> 17) == (end >> 17): return ((1 << 12) - 1) // 7 + (beg >> 17)
    if (beg >> 20) == (end >> 20): return ((1 << 9) - 1) // 7 + (beg >> 20)
    if (beg >> 23) == (end >> 23): return ((1 << 6) - 1) // 7 + (beg >> 23)
    if (beg >> 26) == (end >> 26): return ((1 << 3) - 1) // 7 + (beg >> 26)
    return 0


def _write_bai(bam_path: Path, ref_lengths, alignments):
    """Escribe un índice .bai mínimo pero válido."""
    bai_path = Path(str(bam_path) + ".bai")
    # Agrupar por tid
    by_tid = {i: [] for i in range(len(ref_lengths))}
    n_no_coor = 0
    for a in alignments:
        if a["tid"] < 0:
            n_no_coor += 1
            continue
        by_tid[a["tid"]].append(a)

    with open(bai_path, "wb") as f:
        f.write(b"BAI\x01")
        f.write(struct.pack("<i", len(ref_lengths)))
        for tid in range(len(ref_lengths)):
            alns = by_tid[tid]
            if not alns:
                f.write(struct.pack("<i", 0))   # n_bins
                f.write(struct.pack("<i", 0))   # n_intv
                continue
            # Bins
            bin_map = {}  # bin → list of (chunk_beg, chunk_end)
            linear_idx = {}  # window → min vo
            ref_len = ref_lengths[tid][1]
            n_windows = (ref_len + 16383) // 16384
            for a in alns:
                pos = a["pos"]
                end = pos + max(1, a["ref_len"])
                b = _reg2bin(pos, end)
                bin_map.setdefault(b, []).append((a["_vo_start"], a["_vo_end"]))
                w_start = pos >> 14
                w_end = (end - 1) >> 14
                for w in range(w_start, w_end + 1):
                    if w not in linear_idx or a["_vo_start"] < linear_idx[w]:
                        linear_idx[w] = a["_vo_start"]
            # Coalescer chunks contiguos por bin
            f.write(struct.pack("<i", len(bin_map)))
            for b, chunks in bin_map.items():
                chunks.sort()
                # merge
                merged = [chunks[0]]
                for cb, ce in chunks[1:]:
                    lb, le = merged[-1]
                    if cb <= le:
                        merged[-1] = (lb, max(le, ce))
                    else:
                        merged.append((cb, ce))
                f.write(struct.pack("<I", b))
                f.write(struct.pack("<i", len(merged)))
                for cb, ce in merged:
                    f.write(struct.pack("<QQ", cb, ce))
            # Linear index
            f.write(struct.pack("<i", n_windows))
            last_vo = 0
            for w in range(n_windows):
                vo = linear_idx.get(w, last_vo)
                last_vo = vo
                f.write(struct.pack("<Q", vo))
        # n_no_coor (opcional)
        f.write(struct.pack("<Q", n_no_coor))
