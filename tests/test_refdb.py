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
from fakedb import FakeConn, make_bam  # noqa: E402

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

    def test_existing_read_set_keeps_its_runs(self):
        # OG838 (migration 093): the 240813 read set is the 240813 cell alone; the date rule also
        # gives the 240520 cell, and a GenomeScope re-push must not add it back.
        conn = FakeConn(pacbio_rows={"OG838": [("OG838_A_PACB_240520", "240520"), ("OG838_B_PACB_240813", "240813")]},
                        existing_read_sets={("OG838", "240813"): (206, True)})
        refdb.upsert_read_set(conn.cursor(), "ref_genomescope", "OG838", "240813", {"genomesize": 6})
        self.assertFalse(conn.sql_matching("INSERT INTO ref_read_set ("))
        self.assertFalse(conn.sql_matching("INSERT INTO ref_read_set_run"))
        self.assertEqual(conn.inserts("ref_genomescope")[0]["read_set_id"], 206)

    def test_existing_read_set_without_runs_gets_them(self):
        conn = FakeConn(pacbio_rows={"OG16": [("OG16G_D_SBL_PACB_240427_LA", "240427")]},
                        existing_read_sets={("OG16", "240427"): (55, False)})
        refdb.upsert_read_set(conn.cursor(), "ref_genomescope", "OG16", "240417", {"genomesize": 6})
        runs = [p for s, p in conn.statements if s.startswith("INSERT INTO ref_read_set_run")]
        self.assertEqual(runs, [(101, "OG16G_D_SBL_PACB_240427_LA")])

    def test_delete_keeps_a_read_set_an_assembly_uses(self):
        conn = FakeConn(pacbio_rows={"OG762": [("OG762H_D_SBL_PACB_241018_LA", "241018")]})
        refdb.upsert_read_set(conn.cursor(), "ref_genomescope", "OG762", "241019", {"genomesize": None})
        self.assertIn("NOT EXISTS (SELECT 1 FROM ref_assembly a WHERE a.read_set_id = rs.read_set_id)",
                      conn.sql_matching("DELETE FROM ref_read_set rs")[0])


CELLS = [("OG90G_D_SBL_PACB_240821_LA", "240821", "bc2008"),
         ("OG90G_D_SBL2_PACB_240821_LA", "240821", "bc2009"),
         ("OG90G_D_SBL_PACB_250301_LA", "250301", "bc2008")]


class AdapterFilterRuns(unittest.TestCase):
    def test_parse_stats_and_names(self):
        self.assertEqual(refdb.parse_adapter_stats("Number of adapter contaminated ccs reads: 12 (0.01%)"), 12)
        self.assertIsNone(refdb.parse_adapter_stats("nothing here"))
        self.assertEqual(refdb.adapter_stats_cell("x/OG1_m84154_241127_054344_s3.hifi_reads.bc2025.stats"),
                         ("m84154_241127_054344_s3", "241127", "bc2025"))
        self.assertEqual(refdb.adapter_stats_cell("m64497e_230525_065409.hifi_reads.bc2009--bc2009.stats"),
                         ("m64497e_230525_065409", "230525", "bc2009"))
        self.assertEqual(refdb.adapter_stats_cell("m64497e_230525_065409.unbarcoded.hifi_reads.stats")[2], "none")
        self.assertEqual(refdb.adapter_stats_cell("SSB_cell1_ccs.stats"), (None, None, "none"))

    def test_barcode_then_nearest_date_and_unassigned_follows(self):
        files = [("a/OG90_m84154_240822_054344_s1.hifi_reads.bc2009.stats", 4),
                 ("a/OG90_m84154_240822_054344_s1.hifi_reads.unassigned.stats", 1),
                 ("b/OG90_m84154_240822_054344_s1.hifi_reads.bc2009.stats", 4),        # a copy: once
                 ("a/OG90_m84154_250302_054344_s2.hifi_reads.bc2008.stats", 7)]
        per_run, left_out = refdb.map_adapter_stats(files, CELLS)
        self.assertEqual(left_out, [])
        self.assertEqual({k: v[0] for k, v in per_run.items()},
                         {"OG90G_D_SBL2_PACB_240821_LA": 5, "OG90G_D_SBL_PACB_250301_LA": 7})

    def test_disagreeing_copies_raise(self):
        files = [("a/m84154_240822_054344_s1.hifi_reads.bc2009.stats", 4),
                 ("b/m84154_240822_054344_s1.hifi_reads.bc2009.stats", 5)]
        with self.assertRaises(ValueError):
            refdb.map_adapter_stats(files, CELLS)

    def test_ambiguous_files_left_out(self):
        files = [("m84154_240822_054344_s1.hifi_reads.bc2008.stats", 0),
                 ("m84154_240822_054344_s1.hifi_reads.bc2009.stats", 0),
                 ("m84154_240822_054344_s1.hifi_reads.unassigned.stats", 3),
                 ("m84154_231001_054344_s1.hifi_reads.bc2008.stats", 2),             # no row within 21 days
                 ("SSB_cell1_ccs.stats", 0)]
        per_run, left_out = refdb.map_adapter_stats(files, CELLS)
        self.assertEqual(set(per_run), {"OG90G_D_SBL_PACB_240821_LA", "OG90G_D_SBL2_PACB_240821_LA"})
        self.assertEqual(len(left_out), 3)
        # A zero-count unassigned file of such a movie adds nothing and is dropped quietly.
        files[2] = ("m84154_240822_054344_s1.hifi_reads.unassigned.stats", 0)
        self.assertEqual(len(refdb.map_adapter_stats(files, CELLS)[1]), 2)

    def test_write_adapter_filter_runs_then_read_set_in_one_transaction(self):
        conn = FakeConn(pacbio_cells={"OG90": CELLS},
                        pacbio_rows={"OG90": [("OG90G_D_SBL2_PACB_240821_LA", "240821")]})
        files = {"OG90": [("OG90_m84154_240822_054344_s1.hifi_reads.bc2009.stats", 4)]}
        with mock.patch.object(refdb, "connect", return_value=conn), redirect_stdout(io.StringIO()):
            refdb.write_adapter_filter({}, [("OG90", "240821", {"contam_reads": 4})], files)
        self.assertEqual((conn.commits, conn.rollbacks), (1, 0))
        runs = [p for s, p in conn.statements if s.startswith("INSERT INTO ref_hifi_adapter_filter_run")]
        self.assertEqual(runs, [("OG90G_D_SBL2_PACB_240821_LA", 4, "OG90_m84154_240822_054344_s1.hifi_reads.bc2009.stats")])
        self.assertEqual(conn.inserts("ref_hifi_adapter_filter")[0]["contam_reads"], 4)

    def test_negative_count_rejected(self):
        with self.assertRaises(ValueError):
            refdb.upsert_adapter_filter_run(FakeConn().cursor(), "OG90G_D_SBL2_PACB_240821_LA", -1)

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


class HicRuns(unittest.TestCase):
    """Migration 094: an assembly's Hi-C runs are the flowcells in its Hi-C BAM's read names."""

    def bam(self, names, **kw):
        import os
        import tempfile
        d = tempfile.mkdtemp()
        self.addCleanup(__import__("shutil").rmtree, d)
        path = os.path.join(d, "OG2119_v260408.hic1.hap1.mapped.PT.bam")
        make_bam(path, names, **kw)
        return path

    def test_flowcells_counted_across_blocks(self):
        names = ([f"A01995:73:HGK2LDRX7:1:1101:{i}:1000" for i in range(2500)]
                 + [f"A01995:77:HNTVNDRX7:2:2202:{i}:2000" for i in range(700)])
        self.assertEqual(refdb.bam_flowcells(self.bam(names, block_records=300)),
                         {"HGK2LDRX7": 2500, "HNTVNDRX7": 700})

    def test_names_without_a_flowcell_are_skipped(self):
        path = self.bam(["SRR1234.1", "A01995:73:HGK2LDRX7:1:1101:1:1"])
        self.assertEqual(refdb.bam_flowcells(path), {"HGK2LDRX7": 1})

    def test_only_the_first_bytes_are_read(self):
        names = [f"A01995:73:HGK2LDRX7:1:1101:{i}:1000" for i in range(5000)]
        n = refdb.bam_flowcells(self.bam(names, block_records=100), max_bytes=2000)["HGK2LDRX7"]
        self.assertTrue(0 < n < 5000)

    def test_not_a_bam_or_no_flowcell_raises(self):
        import os
        import tempfile
        d = tempfile.mkdtemp()
        self.addCleanup(__import__("shutil").rmtree, d)
        txt = os.path.join(d, "x.bam")
        with open(txt, "w") as fh:
            fh.write("not a bam")
        with self.assertRaises(ValueError):
            refdb.bam_flowcells(txt)
        with self.assertRaises(ValueError):
            refdb.bam_flowcells(self.bam(["SRR1234.1"]))

    def test_set_hic_flowcells_sorts_and_returns_problems(self):
        conn = FakeConn(assemblies={("OG2119", "260408", "hic1"): 7}, hic_problems={7: ["flowcell X: no run"]})
        out = refdb.set_hic_flowcells(conn.cursor(), "OG2119", "260408", "hic1", {"HNTVNDRX7": 2, "HGK2LDRX7": 5})
        self.assertEqual(out, ["flowcell X: no run"])
        sql, params = conn.statements[-1]
        self.assertTrue(sql.startswith("SELECT set_ref_assembly_hic_flowcells"))
        self.assertEqual(params, (7, ["HGK2LDRX7", "HNTVNDRX7"], [5, 2]))

    def test_set_hic_flowcells_needs_the_assembly(self):
        with self.assertRaises(ValueError):
            refdb.set_hic_flowcells(FakeConn().cursor(), "OG2119", "260408", "hic9", {"HGK2LDRX7": 1})

    def test_write_commits_then_raises_on_problems(self):
        conn = FakeConn(assemblies={("OG1", "240101", "hic1"): 1, ("OG2", "240101", "hic1"): 2},
                        hic_problems={2: ["flowcell HWLY7DRX5 (NOVA_250505_AD): no non-failed Hi-C sequencing row of OG2"]})
        items = [(("OG1", "240101", "hic1"), {"HGK2LDRX7": 1}), (("OG2", "240101", "hic1"), {"HWLY7DRX5": 1})]
        with mock.patch.object(refdb, "connect", return_value=conn), redirect_stdout(io.StringIO()) as out, \
                self.assertRaises(RuntimeError) as cm:
            refdb.write_hic_flowcells({}, items)
        self.assertEqual((conn.commits, conn.rollbacks), (1, 0))
        self.assertIn("1 linked, 1 unresolved", out.getvalue())
        self.assertIn("OG2 240101.hic1: flowcell HWLY7DRX5", str(cm.exception))


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
