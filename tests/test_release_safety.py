from __future__ import annotations

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {".py", ".r", ".md", ".txt", ".cff", ".json", ".yml", ".yaml", ".example"}
FORBIDDEN_SUFFIXES = {".csv", ".tsv", ".dta", ".sav", ".sas7bdat", ".xpt", ".rds", ".rdata", ".doc", ".docx", ".pdf", ".xlsx"}
SENSITIVE_PATTERNS = [
    re.compile(r"[A-Za-z]:[\\/](?:Users|WEST|Codex-outcoms)", re.I),
    re.compile(r"(?:ghp|gho|github_pat)_[A-Za-z0-9_]+"),
    re.compile(r"\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b", re.I),
    re.compile(r"(?<!\d)(?:\+?86[ -]?)?1[3-9]\d{9}(?!\d)"),
]


class ReleaseSafetyTest(unittest.TestCase):
    def test_no_restricted_file_types(self) -> None:
        bad = [str(path.relative_to(ROOT)) for path in ROOT.rglob("*") if path.is_file() and path.suffix.lower() in FORBIDDEN_SUFFIXES]
        self.assertEqual(bad, [])

    def test_no_sensitive_literals(self) -> None:
        findings: list[str] = []
        for path in ROOT.rglob("*"):
            if not path.is_file() or ".git" in path.parts or path.suffix.lower() not in TEXT_SUFFIXES:
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
            for pattern in SENSITIVE_PATTERNS:
                if pattern.search(text):
                    findings.append(f"{path.relative_to(ROOT)}: {pattern.pattern}")
        self.assertEqual(findings, [])


if __name__ == "__main__":
    unittest.main()
