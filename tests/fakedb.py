"""A stand-in for a psycopg2 connection: records SQL instead of running it (Phase RG tests).

Shared, identical, between OceanOmics-OceanGenomes-ref-genomes and -post-curation.
"""
import re


class FakeCursor:
    def __init__(self, conn):
        self.conn = conn
        self.rowcount = 0
        self._result = None
        self._rows = []

    def execute(self, sql, params=None):
        flat = " ".join(sql.split())
        self.conn.statements.append((flat, params))
        if self.conn.fail_on and self.conn.fail_on in flat:
            raise RuntimeError(f"simulated database error on: {self.conn.fail_on}")
        up = flat.upper()
        if "RETURNING" in up:
            self.conn.next_id += 1
            self._result = (self.conn.next_id,)
        elif up.startswith("SELECT P.PART_ID"):
            self._result = self.conn.existing_part
        elif up.startswith("SELECT S.SEQUENCING_ID, S.SEQ_DATE FROM V_SEQUENCING_BROWSE"):
            self._rows = list(self.conn.pacbio_rows.get(params[0], []))
        elif up.startswith("SELECT DISTINCT S.SEQUENCING_ID, S.SEQ_DATE, Q.BARCODE"):
            self._rows = list(self.conn.pacbio_cells.get(params[0], []))
        elif up.startswith("SELECT RS.READ_SET_ID"):
            self._result = self.conn.existing_read_sets.get(params)
        elif up.startswith("SELECT ASSEMBLY_ID FROM REF_ASSEMBLY WHERE"):
            aid = self.conn.assemblies.get(tuple(params))
            self._result = (aid,) if aid else None
        elif up.startswith("SELECT SET_REF_ASSEMBLY_HIC_FLOWCELLS"):
            self._result = (list(self.conn.hic_problems.get(params[0], [])),)
        elif up.startswith("DELETE"):
            self.rowcount = self.conn.delete_rowcount
        else:
            self._result = None

    def fetchone(self):
        return self._result

    def fetchall(self):
        return self._rows

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class FakeConn:
    def __init__(self, existing_part=None, delete_rowcount=1, fail_on=None, pacbio_rows=None,
                 pacbio_cells=None, existing_read_sets=None, assemblies=None, hic_problems=None):
        self.statements = []
        self.commits = 0
        self.rollbacks = 0
        self.closed = False
        self.next_id = 100
        self.existing_part = existing_part      # (part_id, assembly_id) returned by find_part
        self.delete_rowcount = delete_rowcount
        self.fail_on = fail_on                  # substring of SQL that raises
        self.pacbio_rows = pacbio_rows or {}    # og_id -> [(sequencing_id, seq_date)] for read sets
        self.pacbio_cells = pacbio_cells or {}  # og_id -> [(sequencing_id, seq_date, barcode)]
        # (og_id, seq_date) -> (read_set_id, has runs) for read sets already in the database
        self.existing_read_sets = existing_read_sets or {}
        self.assemblies = assemblies or {}      # (og_id, seq_date, version) -> assembly_id
        self.hic_problems = hic_problems or {}  # assembly_id -> what set_ref_assembly_hic_flowcells reports

    def cursor(self):
        return FakeCursor(self)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if exc_type:
            self.rollbacks += 1
        else:
            self.commits += 1
        return False

    def close(self):
        self.closed = True

    # -- helpers for assertions --
    def inserts(self, table):
        """The parameter dicts of every INSERT INTO <table>."""
        return [p for s, p in self.statements if re.match(rf"INSERT INTO {table} \(", s)]

    def tables_inserted(self):
        out = {}
        for s, _ in self.statements:
            m = re.match(r"INSERT INTO (\w+) ", s)
            if m:
                out[m.group(1)] = out.get(m.group(1), 0) + 1
        return out

    def sql_matching(self, needle):
        return [s for s, _ in self.statements if needle in s]


def make_bam(path, read_names, block_records=1000):
    """Write a minimal coordinate-free BAM (one reference, unmapped records with only a name), in
    BGZF blocks of block_records records, as the Hi-C BAM tests need."""
    import struct
    import zlib

    def bgzf(raw):
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        body = c.compress(raw) + c.flush()
        header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
        return (header + struct.pack("<H", len(body) + 25) + body
                + struct.pack("<II", zlib.crc32(raw) & 0xffffffff, len(raw)))

    text = b"@HD\tVN:1.5\tSO:coordinate\n"
    head = b"BAM\x01" + struct.pack("<i", len(text)) + text + struct.pack("<i", 1)
    head += struct.pack("<i", 5) + b"chr1\x00" + struct.pack("<i", 1000)
    blocks = [bgzf(head)]
    recs = b""
    for i, name in enumerate(read_names, 1):
        n = name.encode() + b"\x00"
        body = struct.pack("<iiBBHHHiiii", 0, i, len(n), 60, 4680, 0, 77, 0, -1, -1, 0) + n
        recs += struct.pack("<i", len(body)) + body
        if i % block_records == 0:
            blocks.append(bgzf(recs))
            recs = b""
    if recs:
        blocks.append(bgzf(recs))
    blocks.append(bgzf(b""))
    with open(path, "wb") as fh:
        fh.write(b"".join(blocks))
