"""06_decontam_compile.py and 06a_push_decontam_results_to_sqldb.py (genomes_db_review migration 083).

The report shapes are those of OG2089 v260721.hic1 on scratch, 2026-11-09. Run: python -m unittest discover tests
"""
import csv
import io
import runpy
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from test_push_scripts import COMPILE, PushScriptCase

COMPILER = COMPILE / "06_decontam_compile.py"
PUSH = COMPILE / "06a_push_decontam_results_to_sqldb.py"
compiler = runpy.run_path(str(COMPILER))

TIARA = """sequence_id\tclass_fst_stage\tclass_snd_stage\torg\tbac\tarc\teuk\tunk1\tpla\tunk2\tmit
scaffold_1\teukarya\tn/a\t0\t0\t0\t1\t0\tn/a\tn/a\tn/a
scaffold_2\torganelle\tmitochondrion\t1\t0\t0\t0\t0\tn/a\tn/a\t1
scaffold_3\tbacteria\tn/a\t0\t1\t0\t0\t0\tn/a\tn/a\tn/a
scaffold_4\torganelle\tunknown\t1\t0\t0\t0\t0\tn/a\tn/a\tn/a
scaffold_5\tunknown\tn/a\t0\t0\t0\t0\t1\tn/a\tn/a\tn/a
scaffold_6\tarchaea\tn/a\t0\t0\t1\t0\t0\tn/a\tn/a\tn/a
"""
FCS = """##[["FCS genome report", 2, 1], {}]
#seq_id\tstart_pos\tend_pos\tseq_len\taction\tdiv\tagg_cont_cov\ttop_tax_name
scaffold_7\t1\t500\t500\tEXCLUDE\tanml:molluscs\t64\tOstrea edulis
scaffold_1\t1\t100\t5000\tTRIM\tprok:bacteria\t2\tVibrio
scaffold_1\t4901\t5000\t5000\tTRIM\tprok:bacteria\t2\tVibrio
"""
LENGTHS = {"scaffold_1": 5000, "scaffold_2": 16000, "scaffold_3": 2000, "scaffold_4": 1200,
           "scaffold_5": 1500, "scaffold_6": 1100, "scaffold_7": 500}


def fasta(lengths, wrap=60):
    out = []
    for name, n in lengths.items():
        out.append(f">{name} extra words\n")
        seq = "A" * n
        out += [seq[i:i + wrap] + "\n" for i in range(0, n, wrap)]
    return "".join(out)


class Compile(unittest.TestCase):
    def scratch(self, fasta_lengths=LENGTHS, extra=None):
        base = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, base)
        og = base / "OG2089"
        (og / "05-decontamination/NCBI/out").mkdir(parents=True)
        (og / "05-decontamination/tiara").mkdir(parents=True)
        (og / "04-scaffolding/yahs").mkdir(parents=True)
        for hap in ("hap1", "hap2"):
            (og / f"05-decontamination/tiara/OG2089_v260721.hic1.2.yahs.{hap}.tiara.txt").write_text(TIARA)
            (og / f"05-decontamination/tiara/log_OG2089_v260721.hic1.2.yahs.{hap}.tiara.txt").write_text("log")
            (og / f"05-decontamination/NCBI/out/OG2089_v260721.hic1_.2.yahs.{hap}.NCBI.7941.fcs_gx_report.txt").write_text(FCS)
            (og / f"04-scaffolding/yahs/OG2089_v260721.hic1.1.yahs.{hap}_scaffolds_final.fa").write_text(fasta(fasta_lengths))
        for rel, text in (extra or {}).items():
            (base / rel).write_text(text)
        return base

    def run_compile(self, base):
        out = base / "report.tsv"
        res = subprocess.run([sys.executable, str(COMPILER), "--base", str(base), "--out", str(out)],
                             capture_output=True, text=True, check=True)
        return list(csv.DictReader(open(out), delimiter="\t")), res.stderr

    def test_counts_and_bp(self):
        rows, _ = self.run_compile(self.scratch())
        self.assertEqual([(r["og_id"], r["seq_date"], r["version"], r["haplotype"]) for r in rows],
                         [("OG2089", "260721", "hic1", "hap1"), ("OG2089", "260721", "hic1", "hap2")])
        r = rows[0]
        self.assertEqual((r["input_seqs"], r["input_bp"]), ("7", str(sum(LENGTHS.values()))))
        # The module's rule: mitochondrion, plastid, prokarya/bacteria/archaea; not organelle/unknown.
        self.assertEqual((r["tiara_mito_seqs"], r["tiara_mito_bp"]), ("1", "16000"))
        self.assertEqual((r["tiara_plastid_seqs"], r["tiara_plastid_bp"]), ("0", "0"))
        self.assertEqual((r["tiara_prokaryote_seqs"], r["tiara_prokaryote_bp"]), ("2", "3100"))
        self.assertEqual((r["fcs_exclude_seqs"], r["fcs_exclude_bp"]), ("1", "500"))
        self.assertEqual((r["fcs_trim_seqs"], r["fcs_trim_bp"]), ("1", "200"))
        self.assertEqual((r["fcs_fix_seqs"], r["fcs_review_seqs"]), ("0", "0"))

    def test_contig_missing_from_fasta_blanks_bp(self):
        short = {k: v for k, v in LENGTHS.items() if k != "scaffold_3"}
        rows, err = self.run_compile(self.scratch(fasta_lengths=short))
        self.assertEqual(rows[0]["tiara_prokaryote_seqs"], "2")
        self.assertEqual(rows[0]["tiara_prokaryote_bp"], "")
        self.assertEqual(rows[0]["tiara_mito_bp"], "")
        self.assertIn("1 Tiara contig(s) not in the scaffold FASTA", err)

    def test_post_tiara_fasta_decides_what_was_removed(self):
        # Older runs flagged contigs without removing them (OG824 v241004.hic1 on Acacia).
        row, problem = compiler["compile_row"](
            ("OG2089", "260721", "hic1", "hap1"), TIARA.splitlines(True), FCS.splitlines(True),
            fasta(LENGTHS).splitlines(True),
            fasta({k: v for k, v in LENGTHS.items() if k != "scaffold_2"}).splitlines(True))
        self.assertIsNone(problem)
        self.assertEqual((row["tiara_mito_seqs"], row["tiara_mito_bp"]), (1, 16000))
        self.assertEqual((row["tiara_prokaryote_seqs"], row["tiara_prokaryote_bp"]), (0, 0))

    def test_index_prefers_stage_2_and_final_scaffolds(self):
        files = compiler["index"]([
            "OG90/OG90_v240821.hic1/decontamination/tiara/OG90_v240821.hic1.1.yahs.hap1.tiara.txt",
            "OG90/OG90_v240821.hic1/decontamination/tiara/OG90_v240821.hic1.2.yahs.hap1.tiara.txt",
            "OG90/OG90_v240821.hic1/assembly/OG90_v240821.hic1.2.yahs.hap1.fa",
            "OG90/OG90_v240821.hic1/assembly/OG90_v240821.hic1.2.yahs.hap1_scaffolds_final.fa",
            "OG90/OG90_v240821.hic1/assembly/OG90_v240821.hic1.2.tiara.hap1_scaffolds.fa",
            "OG860/OG860_240927.hic1/decontamination/tiara/OG860_240927.hic1.2.yahs.hap1.tiara.txt",
            "OG90/OG90_bc2008_v240821.hic1/decontamination/tiara/OG90_bc2008_v240821.hic1.2.yahs.hap1.tiara.txt",
        ])
        f = files[("OG90", "240821", "hic1", "hap1")]
        self.assertTrue(f["tiara"].endswith(".2.yahs.hap1.tiara.txt"))
        self.assertTrue(f["fasta"].endswith("_scaffolds_final.fa"))
        self.assertIn(("OG860", "240927", "hic1", "hap1"), files)        # date without the v
        self.assertTrue(f["tiara_out"].endswith(".2.tiara.hap1_scaffolds.fa"))
        self.assertEqual({k[0] for k in files}, {"OG90", "OG860"})       # OG90_bc2008 is not an OG


class Push(PushScriptCase):
    def test_writes_stage_2_in_one_transaction(self):
        tsv = ("og_id\tseq_date\tversion\thaplotype\tinput_seqs\tinput_bp\tfcs_exclude_seqs\tfcs_exclude_bp\t"
               "fcs_trim_seqs\tfcs_trim_bp\tfcs_fix_seqs\tfcs_review_seqs\ttiara_mito_seqs\ttiara_mito_bp\t"
               "tiara_plastid_seqs\ttiara_plastid_bp\ttiara_prokaryote_seqs\ttiara_prokaryote_bp\n"
               "OG2089\t260721\thic1\thap1\t7\t27300\t1\t500\t1\t200\t0\t0\t1\t16000\t0\t0\t2\t3100\n"
               "OG2089\t260721\thic1\thap2\t7\t27300\t0\t0\t0\t0\t0\t0\t1\t16000\t0\t0\t0\t\n")
        conn, code, out = self.run_script(PUSH, [], ["{tmp}/final_decontam_report.tsv"],
                                          tsv_override={"final_decontam_report.tsv": tsv})
        self.assert_one_commit(conn, code, out)
        parts = [p for s, p in conn.statements if s.startswith("INSERT INTO ref_assembly_part")]
        self.assertEqual({(p[1], p[2]) for p in parts}, {(2, "hap1"), (2, "hap2")})
        rows = conn.inserts("ref_decontamination")
        self.assertEqual([r["tiara_prokaryote_bp"] for r in rows], [3100, None])
        self.assertEqual(rows[0]["fcs_trim_bp"], 200)


if __name__ == "__main__":
    unittest.main()
