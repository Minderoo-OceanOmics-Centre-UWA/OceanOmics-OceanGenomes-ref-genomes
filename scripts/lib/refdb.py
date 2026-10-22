"""refdb: database writes for the reference-genome push scripts (Phase RG, schema v4).

The same file is vendored in OceanOmics-OceanGenomes-ref-genomes and
OceanOmics-OceanGenomes-post-curation as scripts/lib/refdb.py. Change both together:
tests/test_refdb_copies.py fails when they differ.

The tables are the relational design in genomes_db_review/docs/ref_genomes_runbook_rg0_rg5.md:
one row per assembly (ref_assembly), per part of it (ref_assembly_part: stage x haplotype), and
one table per tool keyed on the part. A tool that did not run has no row, not a row of NULLs.

Table names are unqualified, so search_path decides where they resolve: v4 during the Phase RG
window (PGOPTIONS='-c search_path=v4,public'), public after the RG.6 swap.

Every write runs in one transaction per input file, and any error exits non-zero, so a caller
such as post_pipeline.sh or post-hichifi-assembly.py sees the failure.
"""

import configparser
import math
import re
import sys
from pathlib import Path

# Columns each tool table accepts. Anything else is a programming error.
TOOL_TABLES = {
    "ref_gfastats": ("num_contigs", "contig_n50", "contig_n50_size_mb", "num_scaffolds", "scaffold_n50",
                     "scaffold_n50_size_mb", "largest_scaffold", "largest_scaffold_size_mb",
                     "total_scaffold_length", "total_scaffold_length_size_mb", "gc_content_percent",
                     "num_gaps"),
    "ref_busco": ("dataset", "complete", "single_copy", "multi_copy", "fragmented", "missing", "n_markers",
                  "internal_stop_codon_percent", "scaffold_n50_bus", "contigs_n50_bus", "percent_gaps",
                  "number_of_scaffolds"),
    "ref_merqury_qv": ("unique_k_mers_assembly", "k_mers_total", "qv", "error"),
    "ref_merqury_completeness": ("solid_k_mers", "total_k_mers", "completeness"),
    "ref_omnic": ("total", "total_unmapped", "total_single_sided_mapped", "total_mapped", "total_dups",
                  "total_nodups", "cis", "trans"),
    "ref_seqkit_stats": ("format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len"),
    "ref_chromosome_stats": ("num_chromosomes", "pct_assigned", "pct_no_super", "num_seq_no_super",
                             "max_len_no_super"),
}
READ_SET_TABLES = {
    "ref_genomescope": ("homozygosity", "heterozygosity", "genomesize", "repeatsize", "uniquesize",
                        "modelfit", "errorrate"),
    "ref_hifi_adapter_filter": ("contam_reads",),
}
HIFI_READ_QC_KEY = ("og_id", "tissue", "ext_type", "lib_code", "run_id")
HIFI_READ_QC_COLS = ("barcode", "barcode_quality", "hifi_reads", "hifi_read_length", "hifi_read_quality",
                     "hifi_yield", "polymerase_read_length")
HAPLOTYPES = ("hap1", "hap2", "dual", "p_ctg", "a_ctg")
SEQ_DATE_RE = re.compile(r"^[0-9]{6}$")
VERSION_RE = re.compile(r"^(hic|hifi|ont)[0-9]+$")


# -- input parsing --------------------------------------------------------------------------------

def load_db_config(config_file):
    """psycopg2.connect() keyword arguments from the [postgres] section of an .ini file."""
    path = Path(config_file).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"Config file '{path}' does not exist.")
    config = configparser.ConfigParser()
    config.read(path)
    if not config.has_section("postgres"):
        raise ValueError("Missing [postgres] section in config file.")
    missing = [k for k in ("dbname", "user", "password", "host", "port") if not config.has_option("postgres", k)]
    if missing:
        raise ValueError(f"Missing {missing} in [postgres] section of config file.")
    return {
        "dbname": config.get("postgres", "dbname"),
        "user": config.get("postgres", "user"),
        "password": config.get("postgres", "password"),
        "host": config.get("postgres", "host"),
        "port": config.getint("postgres", "port"),
    }


def _blank(val):
    if val is None:
        return True
    if isinstance(val, float) and math.isnan(val):
        return True
    s = str(val).strip()
    return s == "" or s.lower() == "nan"


def parse_int(val):
    """'1,234' / 1234.0 / '1234' -> 1234; blank, NaN or unparseable -> None."""
    if _blank(val):
        return None
    try:
        return int(float(str(val).strip().replace(",", "")))
    except ValueError:
        return None


def parse_float(val):
    """'12.5%' / '1,234.5' -> float; blank, NaN or unparseable -> None."""
    if _blank(val):
        return None
    try:
        return float(str(val).strip().replace(",", "").replace("%", ""))
    except ValueError:
        return None


def clean(val):
    """A value psycopg2 can bind: numpy scalars to Python, NaN and pandas NA to None."""
    if val is None:
        return None
    if hasattr(val, "item") and not isinstance(val, (str, bytes)):
        val = val.item()
    if isinstance(val, float) and math.isnan(val):
        return None
    if type(val).__name__ in ("NAType", "NaTType"):
        return None
    return val


def load_samplesheet(path, column):
    """{sample: value} for one samplesheet column, skipping blanks. None path -> {}."""
    import csv
    if not path:
        return {}
    out = {}
    with open(Path(path).expanduser(), newline="") as f:
        for row in csv.DictReader(f):
            og = (row.get("sample") or "").strip()
            val = (row.get(column) or "").strip()
            if og and val:
                out[og] = val
    return out


# -- assemblies and parts -------------------------------------------------------------------------

def check_part_key(og_id, seq_date, version, stage, haplotype):
    """Raise ValueError naming the bad part of a key, before the database sees it."""
    problems = []
    if not og_id or not str(og_id).startswith("OG"):
        problems.append(f"og_id {og_id!r}")
    if seq_date is None or not SEQ_DATE_RE.match(str(seq_date)):
        problems.append(f"seq_date {seq_date!r} (want YYMMDD)")
    if version is None or not VERSION_RE.match(str(version)):
        problems.append(f"version {version!r} (want hicN, hifiN or ontN)")
    if stage is None or str(stage) not in ("0", "1", "2", "3"):
        problems.append(f"stage {stage!r} (want 0-3)")
    if haplotype not in HAPLOTYPES:
        problems.append(f"haplotype {haplotype!r} (want one of {', '.join(HAPLOTYPES)})")
    if problems:
        raise ValueError(f"bad assembly key for {og_id}: " + "; ".join(problems))


def ensure_assembly(cur, og_id, seq_date, version):
    cur.execute(
        "INSERT INTO ref_assembly (og_id, seq_date, version) VALUES (%s, %s, %s) "
        "ON CONFLICT (og_id, seq_date, version) DO UPDATE SET og_id = EXCLUDED.og_id "
        "RETURNING assembly_id",
        (og_id, str(seq_date), version))
    return cur.fetchone()[0]


def ensure_part(cur, assembly_id, stage, haplotype):
    cur.execute(
        "INSERT INTO ref_assembly_part (assembly_id, stage, haplotype) VALUES (%s, %s, %s) "
        "ON CONFLICT (assembly_id, stage, haplotype) DO UPDATE SET stage = EXCLUDED.stage "
        "RETURNING part_id",
        (assembly_id, int(stage), haplotype))
    return cur.fetchone()[0]


def find_part(cur, og_id, seq_date, version, stage, haplotype):
    cur.execute(
        "SELECT p.part_id, a.assembly_id FROM ref_assembly a "
        "JOIN ref_assembly_part p ON p.assembly_id = a.assembly_id "
        "WHERE a.og_id = %s AND a.seq_date = %s AND a.version = %s AND p.stage = %s AND p.haplotype = %s",
        (og_id, str(seq_date), version, int(stage), haplotype))
    return cur.fetchone()


def prune(cur, part_id, assembly_id):
    """Drop a part with no tool rows left, then an assembly with no parts left: no empty rows."""
    tool_exists = " OR ".join(f"EXISTS (SELECT 1 FROM {t} x WHERE x.part_id = %(p)s)" for t in TOOL_TABLES)
    cur.execute(f"DELETE FROM ref_assembly_part WHERE part_id = %(p)s AND NOT ({tool_exists})", {"p": part_id})
    cur.execute("DELETE FROM ref_assembly a WHERE a.assembly_id = %s "
                "AND NOT EXISTS (SELECT 1 FROM ref_assembly_part p WHERE p.assembly_id = a.assembly_id)",
                (assembly_id,))


def upsert_metrics(cur, table, key, values):
    """Write one tool's values for one part.

    key: dict with og_id, seq_date, version, stage, haplotype.
    values: dict of that tool's columns. Columns not given are left as they are.
    If every value is None the tool's row for that part is removed (as the old wide-table upsert
    blanked the block), and an emptied part or assembly is removed with it.
    Returns 'upserted', 'deleted' or 'unchanged'.
    """
    if table not in TOOL_TABLES:
        raise ValueError(f"{table} is not a tool table")
    unknown = set(values) - set(TOOL_TABLES[table])
    if unknown:
        raise ValueError(f"{table} has no column(s) {sorted(unknown)}")
    k = {c: clean(key.get(c)) for c in ("og_id", "seq_date", "version", "stage", "haplotype")}
    check_part_key(**k)
    vals = {c: clean(v) for c, v in values.items()}

    if all(v is None for v in vals.values()):
        found = find_part(cur, **k)
        if not found:
            return "unchanged"
        part_id, assembly_id = found
        cur.execute(f"DELETE FROM {table} WHERE part_id = %s", (part_id,))
        prune(cur, part_id, assembly_id)
        return "deleted"

    assembly_id = ensure_assembly(cur, k["og_id"], k["seq_date"], k["version"])
    part_id = ensure_part(cur, assembly_id, k["stage"], k["haplotype"])
    cols = list(vals)
    quoted = [f'"{c}"' if c == "type" else c for c in cols]
    cur.execute(
        f"INSERT INTO {table} (part_id, {', '.join(quoted)}) "
        f"VALUES (%(part_id)s, {', '.join('%(' + c + ')s' for c in cols)}) "
        f"ON CONFLICT (part_id) DO UPDATE SET {', '.join(f'{q} = EXCLUDED.{q}' for q in quoted)}, loaded_at = now()",
        dict(vals, part_id=part_id))
    return "upserted"


# -- read sets and read QC ------------------------------------------------------------------------

# A samplesheet date can be the PacBio movie date, up to ~10 days off the run's seq_date (OG16:
# v240417, run PACB_240427_LA). Migration 062 dated the backfilled read sets with the same rule.
READ_SET_WINDOW_DAYS = 21


def read_set_runs(cur, og_id, sheet_date):
    """The OG's PacBio DNA sequencing rows a samplesheet date covers: every run up to the date plus
    READ_SET_WINDOW_DAYS. Returns (read set seq_date, [sequencing_id]). The read set's seq_date is
    the latest seq_date among them, or the samplesheet date when the OG has none (external reads),
    which v4.v_read_set_run_readiness then lists."""
    cur.execute(
        "SELECT s.sequencing_id, s.seq_date FROM v_sequencing_browse s "
        "WHERE s.og_id = %s AND s.technology IN ('PacBio HIFI', 'PacBio') AND s.library_type = 'pacbio' "
        "AND s.seq_date ~ '^[0-9]{6}$' "
        "AND to_date(s.seq_date, 'YYMMDD') <= to_date(%s, 'YYMMDD') + %s "
        "ORDER BY s.seq_date, s.sequencing_id",
        (og_id, sheet_date, READ_SET_WINDOW_DAYS))
    runs = cur.fetchall() or []
    return max((d for _, d in runs), default=sheet_date), [sid for sid, _ in runs]


def ensure_read_set(cur, og_id, sheet_date):
    """read_set_id for (og_id, samplesheet date), created if missing, with its runs linked."""
    seq_date, runs = read_set_runs(cur, og_id, sheet_date)
    cur.execute(
        "INSERT INTO ref_read_set (og_id, seq_date) VALUES (%s, %s) "
        "ON CONFLICT (og_id, (coalesce(seq_date, ''))) DO UPDATE SET og_id = EXCLUDED.og_id "
        "RETURNING read_set_id",
        (og_id, seq_date))
    read_set_id = cur.fetchone()[0]
    for sid in runs:
        cur.execute("INSERT INTO ref_read_set_run (read_set_id, sequencing_id) VALUES (%s, %s) "
                    "ON CONFLICT DO NOTHING", (read_set_id, sid))
    return read_set_id


def upsert_read_set(cur, table, og_id, seq_date, values):
    """GenomeScope / HiFiAdapterFilt for one HiFi read set (v4.ref_read_set, migration 062).

    seq_date is the samplesheet's date (YYMMDD) and is required: a result whose reads are unknown is
    not written. All-None values remove the tool row, and the read set too once no tool row is left.
    Returns 'upserted', 'deleted' or 'unchanged'.
    """
    if table not in READ_SET_TABLES:
        raise ValueError(f"{table} is not a read-set table")
    unknown = set(values) - set(READ_SET_TABLES[table])
    if unknown:
        raise ValueError(f"{table} has no column(s) {sorted(unknown)}")
    og_id, seq_date = clean(og_id), clean(seq_date)
    if not og_id or not str(og_id).startswith("OG"):
        raise ValueError(f"bad og_id {og_id!r}")
    if seq_date is None:
        raise ValueError(f"no samplesheet date for {og_id}: pass the samplesheet so the read set is known")
    if not SEQ_DATE_RE.match(str(seq_date)):
        raise ValueError(f"bad seq_date {seq_date!r} for {og_id} (want YYMMDD)")
    seq_date = str(seq_date)
    vals = {c: clean(v) for c, v in values.items()}
    if all(v is None for v in vals.values()):
        key, _ = read_set_runs(cur, og_id, seq_date)
        cur.execute(f"DELETE FROM {table} t USING ref_read_set rs "
                    f"WHERE t.read_set_id = rs.read_set_id AND rs.og_id = %s AND rs.seq_date = %s",
                    (og_id, key))
        deleted = cur.rowcount
        cur.execute("DELETE FROM ref_read_set rs WHERE rs.og_id = %s AND rs.seq_date = %s "
                    "AND NOT EXISTS (SELECT 1 FROM ref_genomescope g WHERE g.read_set_id = rs.read_set_id) "
                    "AND NOT EXISTS (SELECT 1 FROM ref_hifi_adapter_filter f WHERE f.read_set_id = rs.read_set_id)",
                    (og_id, key))
        return "deleted" if deleted else "unchanged"
    read_set_id = ensure_read_set(cur, og_id, seq_date)
    cols = list(vals)
    cur.execute(
        f"INSERT INTO {table} (read_set_id, {', '.join(cols)}) "
        f"VALUES (%(read_set_id)s, {', '.join('%(' + c + ')s' for c in cols)}) "
        f"ON CONFLICT (read_set_id) DO UPDATE SET "
        f"{', '.join(f'{c} = EXCLUDED.{c}' for c in cols)}, loaded_at = now()",
        dict(vals, read_set_id=read_set_id))
    return "upserted"


def upsert_hifi_read_qc(cur, row):
    """One PacBio HiFi library-on-run row. Links the one sequencing row of that library on that run,
    rebuilding the tube as og_id||tissue||'_'||ext_type||'_'||lib_code (NULL if none or several)."""
    r = {c: clean(row.get(c)) for c in HIFI_READ_QC_KEY + HIFI_READ_QC_COLS}
    blank = [c for c in HIFI_READ_QC_KEY if not r[c]]
    if blank:
        raise ValueError(f"hifi_read_qc row with blank {blank}: {r}")
    cur.execute(
        f"INSERT INTO hifi_read_qc ({', '.join(HIFI_READ_QC_KEY)}, sequencing_id, {', '.join(HIFI_READ_QC_COLS)}) "
        f"SELECT {', '.join('%(' + c + ')s' for c in HIFI_READ_QC_KEY)}, "
        "(SELECT CASE WHEN count(*) = 1 THEN min(s.sequencing_id) END FROM sequencing s "
        " WHERE s.library_tube_id = %(og_id)s || %(tissue)s || '_' || %(ext_type)s || '_' || %(lib_code)s "
        "   AND s.run_id = %(run_id)s), "
        f"{', '.join('%(' + c + ')s' for c in HIFI_READ_QC_COLS)} "
        f"ON CONFLICT ({', '.join(HIFI_READ_QC_KEY)}) DO UPDATE SET sequencing_id = EXCLUDED.sequencing_id, "
        f"{', '.join(f'{c} = EXCLUDED.{c}' for c in HIFI_READ_QC_COLS)}, loaded_at = now()",
        r)
    return "upserted"


# -- transactions -----------------------------------------------------------------------------------

def connect(db_params):
    import psycopg2
    return psycopg2.connect(**db_params)


def write_all(db_params, writes, label):
    """Run writes(cur) -> {outcome: count} in one transaction. Commit only if it all succeeds."""
    conn = connect(db_params)
    try:
        with conn:
            with conn.cursor() as cur:
                counts = writes(cur)
    finally:
        conn.close()
    summary = ", ".join(f"{n} {k}" for k, n in sorted(counts.items())) or "nothing to write"
    print(f"✅ {label}: {summary} (one transaction, committed)")
    return counts


def write_tool_rows(db_params, table, rows):
    """rows: iterable of (key dict, values dict). One transaction for the whole file."""
    rows = list(rows)

    def writes(cur):
        counts = {}
        for key, values in rows:
            outcome = upsert_metrics(cur, table, key, values)
            counts[outcome] = counts.get(outcome, 0) + 1
        return counts
    return write_all(db_params, writes, table)


def write_read_set_rows(db_params, table, rows):
    """rows: iterable of (og_id, samplesheet date YYMMDD, values dict). One transaction; a row with
    no date fails the whole file."""
    rows = list(rows)

    def writes(cur):
        counts = {}
        for og_id, seq_date, values in rows:
            outcome = upsert_read_set(cur, table, og_id, seq_date, values)
            counts[outcome] = counts.get(outcome, 0) + 1
        return counts
    return write_all(db_params, writes, table)


def write_hifi_read_qc_rows(db_params, rows):
    rows = list(rows)

    def writes(cur):
        for row in rows:
            upsert_hifi_read_qc(cur, row)
        return {"upserted": len(rows)}
    return write_all(db_params, writes, "hifi_read_qc")


def main_guard(main):
    """Run main(); on any error print it and exit 1, so pipelines see the failure."""
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:  # noqa: BLE001 - report everything, then fail
        print(f"❌ Error: {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(1)
