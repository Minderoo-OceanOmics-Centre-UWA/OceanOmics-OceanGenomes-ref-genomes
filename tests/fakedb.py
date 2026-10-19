"""A stand-in for a psycopg2 connection: records SQL instead of running it (Phase RG tests).

Shared, identical, between OceanOmics-OceanGenomes-ref-genomes and -post-curation.
"""
import re


class FakeCursor:
    def __init__(self, conn):
        self.conn = conn
        self.rowcount = 0
        self._result = None

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
        elif up.startswith("DELETE"):
            self.rowcount = self.conn.delete_rowcount
        else:
            self._result = None

    def fetchone(self):
        return self._result

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class FakeConn:
    def __init__(self, existing_part=None, delete_rowcount=1, fail_on=None):
        self.statements = []
        self.commits = 0
        self.rollbacks = 0
        self.closed = False
        self.next_id = 100
        self.existing_part = existing_part      # (part_id, assembly_id) returned by find_part
        self.delete_rowcount = delete_rowcount
        self.fail_on = fail_on                  # substring of SQL that raises

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
