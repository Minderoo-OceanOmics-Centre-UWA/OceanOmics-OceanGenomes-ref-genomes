"""Unit tests for scripts/lib/refdb.py (Phase RG). Run: python -m unittest discover tests

Shared, identical, between OceanOmics-OceanGenomes-ref-genomes and -post-curation.
"""
import io
import sys
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts" / "lib"))
sys.path.insert(0, str(HERE))
import refdb  # noqa: E402
from fakedb import FakeConn  # noqa: E402

KEY = {"og_id": "OG762", "seq_date": "241019", "version": "hic1", "stage": 3, "haplotype": "hap1"}


class Parsing(unittest.TestCase):
    def test_parse_int(self):
        self.assertEqual(refdb.parse_int("1,234"), 1234)
        self.assertEqual(refdb.parse_int(1234.0), 1234)
        self.assertIsNone(refdb.parse_int(""))
        self.assertIsNone(refdb.parse_int("nan"))
        self.assertIsNone(refdb.parse_int(float("nan")))
        self.assertIsNone(refdb.parse_int("abc"))

    def test_parse_float(self):
        self.assertEqual(refdb.parse_float("99.3762%"), 99.3762)
        self.assertEqual(refdb.parse_float("823,456,968"), 823456968.0)
        self.assertIsNone(refdb.parse_float(None))

    def test_clean_numpy_and_nan(self):
        import numpy as np
        self.assertIs(type(refdb.clean(np.int64(5))), int)
        self.assertIsNone(refdb.clean(np.float64("nan")))
        self.assertIsNone(refdb.clean(float("nan")))
        self.assertEqual(refdb.clean("x"), "x")

    def test_check_part_key_names_every_problem(self):
        with self.assertRaises(ValueError) as cm:
            refdb.check_part_key("OG762", "bc2008", "hi1", 7, None)
        msg = str(cm.exception)
        for part in ("seq_date", "version", "stage", "haplotype"):
            self.assertIn(part, msg)

    def test_check_part_key_accepts_int_seq_date(self):
        # pandas reads a separate seq_date column as an integer
        refdb.check_part_key("OG762", 241019, "hic1", "3", "dual")


class UpsertMetrics(unittest.TestCase):
    def test_values_create_assembly_part_and_tool_row(self):
        conn = FakeConn()
        cur = conn.cursor()
        out = refdb.upsert_metrics(cur, "ref_busco", KEY, {"complete": 99.2, "single_copy": 98.9})
        self.assertEqual(out, "upserted")
        self.assertEqual(list(conn.tables_inserted()), ["ref_assembly", "ref_assembly_part", "ref_busco"])
        busco_sql = conn.sql_matching("INSERT INTO ref_busco")[0]
        self.assertIn("ON CONFLICT (part_id) DO UPDATE SET complete = EXCLUDED.complete", busco_sql)
        self.assertNotIn("dataset", busco_sql)  # columns not given are left alone
        self.assertEqual(conn.inserts("ref_busco")[0]["part_id"], 102)

    def test_type_column_is_quoted(self):
        conn = FakeConn()
        refdb.upsert_metrics(conn.cursor(), "ref_seqkit_stats", KEY, {"format": "FASTA", "type": "DNA"})
        self.assertIn('"type" = EXCLUDED."type"', conn.sql_matching("INSERT INTO ref_seqkit_stats")[0])

    def test_all_none_with_existing_part_deletes_and_prunes(self):
        conn = FakeConn(existing_part=(7, 3))
        out = refdb.upsert_metrics(conn.cursor(), "ref_omnic", KEY, {"total": None, "cis": float("nan")})
        self.assertEqual(out, "deleted")
        self.assertEqual(conn.tables_inserted(), {})
        self.assertTrue(conn.sql_matching("DELETE FROM ref_omnic WHERE part_id"))
        self.assertTrue(conn.sql_matching("DELETE FROM ref_assembly_part WHERE part_id"))
        self.assertTrue(conn.sql_matching("DELETE FROM ref_assembly a WHERE a.assembly_id"))

    def test_all_none_without_part_writes_nothing(self):
        conn = FakeConn(existing_part=None)
        out = refdb.upsert_metrics(conn.cursor(), "ref_omnic", KEY, {"total": None})
        self.assertEqual(out, "unchanged")
        self.assertFalse([s for s, _ in conn.statements if not s.startswith("SELECT")])

    def test_unknown_table_or_column_rejected(self):
        cur = FakeConn().cursor()
        with self.assertRaises(ValueError):
            refdb.upsert_metrics(cur, "ref_genomes", KEY, {"qv": 1})
        with self.assertRaises(ValueError):
            refdb.upsert_metrics(cur, "ref_busco", KEY, {"qv": 1})

    def test_bad_key_rejected_before_any_sql(self):
        conn = FakeConn()
        with self.assertRaises(ValueError):
            refdb.upsert_metrics(conn.cursor(), "ref_gfastats", dict(KEY, haplotype=None), {"num_contigs": 1})
        self.assertEqual(conn.statements, [])


class ReadSets(unittest.TestCase):
    def test_sheet_date_resolves_to_runs_and_their_latest_seq_date(self):
        # OG16: the samplesheet says v240417 (movie date); the run is PACB_240427_LA.
        conn = FakeConn(pacbio_rows={"OG16": [("OG16G_D_SBL_PACB_240427_LA", "240427")]})
        refdb.upsert_read_set(conn.cursor(), "ref_genomescope", "OG16", "240417", {"genomesize": 6})
        rs = [p for s, p in conn.statements if s.startswith("INSERT INTO ref_read_set (")]
        self.assertEqual(rs, [("OG16", "240427")])
        runs = [p for s, p in conn.statements if s.startswith("INSERT INTO ref_read_set_run")]
        self.assertEqual(runs, [(101, "OG16G_D_SBL_PACB_240427_LA")])
        self.assertEqual(conn.inserts("ref_genomescope")[0]["read_set_id"], 101)
        self.assertIn("ON CONFLICT (read_set_id)", conn.sql_matching("INSERT INTO ref_genomescope")[0])

    def test_no_pacbio_rows_keeps_the_sheet_date(self):
        conn = FakeConn()
        refdb.upsert_read_set(conn.cursor(), "ref_hifi_adapter_filter", "OG906", "250212", {"contam_reads": 0})
        rs = [p for s, p in conn.statements if s.startswith("INSERT INTO ref_read_set (")]
        self.assertEqual(rs, [("OG906", "250212")])
        self.assertFalse(conn.sql_matching("INSERT INTO ref_read_set_run"))

    def test_missing_date_rejected_before_any_write(self):
        conn = FakeConn()
        with self.assertRaises(ValueError):
            refdb.upsert_read_set(conn.cursor(), "ref_genomescope", "OG762", None, {"genomesize": 6})
        self.assertEqual(conn.statements, [])

    def test_all_none_deletes_tool_row_then_empty_read_set(self):
        conn = FakeConn(delete_rowcount=1, pacbio_rows={"OG762": [("OG762H_D_SBL_PACB_241018_LA", "241018")]})
        out = refdb.upsert_read_set(conn.cursor(), "ref_hifi_adapter_filter", "OG762", "241019", {"contam_reads": None})
        self.assertEqual(out, "deleted")
        self.assertTrue(conn.sql_matching("DELETE FROM ref_hifi_adapter_filter"))
        self.assertTrue(conn.sql_matching("DELETE FROM ref_read_set rs"))
        self.assertFalse(conn.sql_matching("INSERT INTO"))

    def test_bad_seq_date_rejected(self):
        with self.assertRaises(ValueError):
            refdb.upsert_read_set(FakeConn().cursor(), "ref_genomescope", "OG762", "v241018", {"genomesize": 1})

    def test_hifi_read_qc_resolves_sequencing_and_rejects_blank_key(self):
        conn = FakeConn()
        row = {"og_id": "OG762", "tissue": "H", "ext_type": "D", "lib_code": "SBL2", "run_id": "PACB_241018_LA",
               "hifi_yield": 1}
        refdb.upsert_hifi_read_qc(conn.cursor(), row)
        sql = conn.sql_matching("INSERT INTO hifi_read_qc")[0]
        self.assertIn("count(*) = 1 THEN min(s.sequencing_id)", sql)
        with self.assertRaises(ValueError):
            refdb.upsert_hifi_read_qc(conn.cursor(), dict(row, lib_code=""))


class Transactions(unittest.TestCase):
    def test_one_commit_for_the_whole_file(self):
        conn = FakeConn()
        rows = [(KEY, {"qv": 69.3}), (dict(KEY, haplotype="hap2"), {"qv": 69.2})]
        with mock.patch.object(refdb, "connect", return_value=conn), redirect_stdout(io.StringIO()):
            counts = refdb.write_tool_rows({}, "ref_merqury_qv", rows)
        self.assertEqual(counts, {"upserted": 2})
        self.assertEqual((conn.commits, conn.rollbacks, conn.closed), (1, 0, True))

    def test_error_rolls_back_everything(self):
        conn = FakeConn(fail_on="INSERT INTO ref_merqury_qv")
        with mock.patch.object(refdb, "connect", return_value=conn), self.assertRaises(RuntimeError):
            refdb.write_tool_rows({}, "ref_merqury_qv", [(KEY, {"qv": 1.0})])
        self.assertEqual((conn.commits, conn.rollbacks, conn.closed), (0, 1, True))

    def test_main_guard_exits_1(self):
        def boom():
            raise ValueError("bad row")
        with redirect_stderr(io.StringIO()) as err, self.assertRaises(SystemExit) as cm:
            refdb.main_guard(boom)
        self.assertEqual(cm.exception.code, 1)
        self.assertIn("bad row", err.getvalue())


class Config(unittest.TestCase):
    def test_load_db_config(self):
        import tempfile
        with tempfile.NamedTemporaryFile("w", suffix=".cfg", delete=False) as f:
            f.write("[postgres]\ndbname=d\nuser=u\npassword=p\nhost=h\nport=5432\n")
        self.assertEqual(refdb.load_db_config(f.name)["port"], 5432)
        Path(f.name).write_text("[postgres]\ndbname=d\n")
        with self.assertRaises(ValueError):
            refdb.load_db_config(f.name)
        Path(f.name).unlink()


if __name__ == "__main__":
    unittest.main()
