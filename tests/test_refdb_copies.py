"""scripts/lib/refdb.py and the shared test files must be byte-identical in both pipeline repos.

Shared, identical, between OceanOmics-OceanGenomes-ref-genomes and -post-curation. Skips when the
other repo is not checked out next to this one.
"""
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
SIBLINGS = ["OceanOmics-OceanGenomes-ref-genomes", "OceanOmics-OceanGenomes-post-curation"]
SHARED = ["scripts/lib/refdb.py", "tests/fakedb.py", "tests/test_refdb.py", "tests/test_refdb_copies.py"]


class SharedFilesIdentical(unittest.TestCase):
    def test_identical(self):
        other_name = next(n for n in SIBLINGS if n != REPO.name) if REPO.name in SIBLINGS else None
        other = REPO.parent / other_name if other_name else None
        if not other or not other.is_dir():
            self.skipTest("sibling pipeline repo not checked out next to this one")
        for rel in SHARED:
            with self.subTest(file=rel):
                self.assertEqual((REPO / rel).read_bytes(), (other / rel).read_bytes(),
                                 f"{rel} differs between {REPO.name} and {other.name}; change both")


if __name__ == "__main__":
    unittest.main()
