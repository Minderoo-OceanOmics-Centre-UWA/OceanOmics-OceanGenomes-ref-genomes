"""Phase RG: each push script writes its own v4 table, in one transaction, and fails loudly.

Every script runs as __main__ in a temporary folder holding the compiled TSV it expects, with
refdb.connect replaced by a fake that records SQL (tests/fakedb.py). Fixtures in tests/fixtures are
real compiled results (OceanOmics-Database/SS/ref_genomes_table, OG762 hic1 241019 and others) plus
small hand-made files. Run: python -m unittest discover tests
"""
import csv
import io
import os
import runpy
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
FIX = HERE / "fixtures"
COMPILE = REPO / "scripts" / "03_compile_results"
sys.path.insert(0, str(REPO / "scripts" / "lib"))
sys.path.insert(0, str(HERE))
import refdb  # noqa: E402
from fakedb import FakeConn  # noqa: E402

SHEET_OGS = {r["sample"] for r in csv.DictReader(open(FIX / "samplesheet.csv"))}


def fixture_rows(name, col=0):
    """Data rows of a fixture whose first token (before '_') is an OG in the samplesheet."""
    with open(FIX / name) as f:
        next(f)
        return [line for line in f if line.split("\t")[col].split("_")[0] in SHEET_OGS]


class PushScriptCase(unittest.TestCase):
    def run_script(self, script, files, args=(), conn=None, tsv_override=None):
        conn = conn or FakeConn()
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        for name in files:
            shutil.copy(FIX / name, tmp / name)
        for name, text in (tsv_override or {}).items():
            (tmp / name).write_text(text)
        cfg = tmp / "db.cfg"
        cfg.write_text("[postgres]\ndbname=d\nuser=u\npassword=p\nhost=h\nport=5432\n")
        argv = [str(script), str(cfg)] + [str(a).replace("{fx}", str(FIX)).replace("{tmp}", str(tmp)) for a in args]
        out, cwd, code = io.StringIO(), os.getcwd(), 0
        os.chdir(tmp)
        try:
            with mock.patch.object(sys, "argv", argv), mock.patch.object(refdb, "connect", return_value=conn), \
                    redirect_stdout(out), redirect_stderr(out):
                runpy.run_path(str(script), run_name="__main__")
        except SystemExit as e:
            code = e.code if isinstance(e.code, int) else 1
        finally:
            os.chdir(cwd)
        return conn, code, out.getvalue()

    def assert_one_commit(self, conn, code, out):
        self.assertEqual(code, 0, out[-2000:])
        self.assertEqual((conn.commits, conn.rollbacks), (1, 0))
        self.assertNotIn("INSERT INTO ref_genomes ", "\n".join(s for s, _ in conn.statements))


class ToolTables(PushScriptCase):
    def test_busco_writes_ref_busco_only(self):
        conn, code, out = self.run_script(COMPILE / "03a_push_busco_results_to_sqldb.py",
                                          ["BUSCO_compiled_results.tsv"], ["{fx}/samplesheet.csv"])
        self.assert_one_commit(conn, code, out)
        self.assertEqual(set(conn.tables_inserted()), {"ref_assembly", "ref_assembly_part", "ref_busco"})
        self.assertEqual(len(conn.inserts("ref_busco")), len(fixture_rows("BUSCO_compiled_results.tsv")))
        self.assertIn("Skipping OG85", out)  # in the TSV, not in the samplesheet

    def test_merqury_qv_includes_dual(self):
        conn, code, out = self.run_script(COMPILE / "02b_push_merqury_qv_results_to_sqldb.py",
                                          ["merqury.qv.stats.tsv"], ["{fx}/samplesheet.csv"])
        self.assert_one_commit(conn, code, out)
        self.assertEqual(len(conn.inserts("ref_merqury_qv")), len(fixture_rows("merqury.qv.stats.tsv")))
        haps = {p[2] for p in (pp for s, pp in conn.statements if s.startswith("INSERT INTO ref_assembly_part"))}
        self.assertEqual(haps, {"hap1", "hap2", "dual"})

    def test_merqury_completeness(self):
        conn, code, out = self.run_script(COMPILE / "02d_push_merqury_completeness_results_to_sqldb.py",
                                          ["merqury.completeness.stats.tsv"], ["{fx}/samplesheet.csv"])
        self.assert_one_commit(conn, code, out)
        self.assertEqual([p["completeness"] for p in conn.inserts("ref_merqury_completeness")],
                         [85.2146, 85.1683, 96.0412])

    def test_gfastats(self):
        conn, code, out = self.run_script(COMPILE / "04a_push_gfa_results_to_sqldb.py",
                                          ["final_gfastats_report.txt"], ["{fx}/samplesheet.csv"])
        self.assert_one_commit(conn, code, out)
        self.assertEqual(len(conn.inserts("ref_gfastats")), len(fixture_rows("final_gfastats_report.txt")))

    def test_omnic_maps_stage_names(self):
        conn, code, out = self.run_script(COMPILE / "05a_push_omnic_results_to_sqldb.py",
                                          ["final_omnic_stats_report.txt"], ["{fx}/samplesheet.csv"])
        self.assert_one_commit(conn, code, out)
        parts = [p for s, p in conn.statements if s.startswith("INSERT INTO ref_assembly_part")]
        self.assertEqual({p[1] for p in parts}, {0, 2})
        self.assertEqual(len(conn.inserts("ref_omnic")), len(fixture_rows("final_omnic_stats_report.txt")))

    def test_hifi_gfastats_skips_excluded_variant(self):
        conn, code, out = self.run_script(COMPILE / "04b_hifi_push_gfa_results_to_sqldb.py",
                                          ["hifi_gfastats_report.txt"])
        self.assert_one_commit(conn, code, out)
        self.assertIn("Skipping OG90_bc2008", out)
        asm = [p for s, p in conn.statements if s.startswith("INSERT INTO ref_assembly ")]
        self.assertEqual({a for a in asm}, {("OG762", "241018", "hifi1")})
        self.assertEqual(len(conn.inserts("ref_gfastats")), 2)


class ReadLevel(PushScriptCase):
    def test_genomescope_dated_from_samplesheet(self):
        conn, code, out = self.run_script(COMPILE / "01a_push_genomescope_comp_results_to_sqldb.py",
                                          ["genomescope_compiled_results.tsv"], ["{fx}/samplesheet.csv"])
        self.assert_one_commit(conn, code, out)
        rows = {p["og_id"]: p for p in conn.inserts("ref_genomescope")}
        self.assertEqual(rows["OG762"]["seq_date"], "241019")
        self.assertEqual(rows["OG113"]["homozygosity"], 99.3762)
        self.assertEqual(rows["OG113"]["genomesize"], 823456968)

    def test_adapter_filter_with_and_without_samplesheet(self):
        conn, code, out = self.run_script(COMPILE / "0a_push_hifiadapt_results_to_sqldb.py",
                                          ["hifiadaptorfilt_stats_summary.tsv"], ["{fx}/samplesheet.csv"])
        self.assert_one_commit(conn, code, out)
        self.assertEqual({p["seq_date"] for p in conn.inserts("ref_hifi_adapter_filter")},
                         {"241019", "231028", "240417", "241004"})
        conn, code, out = self.run_script(COMPILE / "0a_push_hifiadapt_results_to_sqldb.py",
                                          ["hifiadaptorfilt_stats_summary.tsv"])
        self.assert_one_commit(conn, code, out)
        self.assertEqual({p["seq_date"] for p in conn.inserts("ref_hifi_adapter_filter")}, {None})

    def test_hifi_yield_sums_smrt_cells(self):
        conn, code, out = self.run_script(COMPILE / "0b_push_hifi_yield_results_to_sqldb.py",
                                          [], ["{fx}/PACB_241018_LA_yield.csv"])
        self.assert_one_commit(conn, code, out)
        rows = {p["og_id"]: p for p in conn.inserts("hifi_read_qc")}
        self.assertEqual(set(rows), {"OG762", "OG799"})
        self.assertEqual(rows["OG762"]["hifi_reads"], 3000000)
        self.assertEqual(rows["OG762"]["barcode"], "bc2032")


class Failures(PushScriptCase):
    def test_database_error_exits_1_and_rolls_back(self):
        conn, code, out = self.run_script(COMPILE / "03a_push_busco_results_to_sqldb.py",
                                          ["BUSCO_compiled_results.tsv"], ["{fx}/samplesheet.csv"],
                                          conn=FakeConn(fail_on="INSERT INTO ref_busco"))
        self.assertEqual(code, 1)
        self.assertEqual((conn.commits, conn.rollbacks), (0, 1))
        self.assertIn("simulated database error", out)

    def test_bad_assembly_name_exits_1_and_commits_nothing(self):
        bad = (FIX / "merqury.qv.stats.tsv").read_text() + "OG762_v241019.hic1.0.hifiasm\t1\t2\t3\t4\n"
        conn, code, out = self.run_script(COMPILE / "02b_push_merqury_qv_results_to_sqldb.py", [],
                                          ["{fx}/samplesheet.csv"], tsv_override={"merqury.qv.stats.tsv": bad})
        self.assertEqual(code, 1)
        self.assertEqual(conn.commits, 0)
        self.assertIn("bad assembly key for OG762: haplotype", out)


if __name__ == "__main__":
    unittest.main()
