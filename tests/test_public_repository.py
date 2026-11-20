"""Checks of this repository as published: the guard passes, the locked files match, and no server path is present."""
import csv
import pathlib
import re
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))


def tracked():
    return subprocess.check_output(["git", "ls-files"], cwd=REPO, text=True).splitlines()


def test_loader_guard_passes():
    import guard

    assert guard.check(REPO)["fixed_commit"] != guard.check(REPO)["pipeline_commit"]


def test_no_server_path_or_account_name():
    needle = re.compile(rb"/" + rb"A\.I_DATA|jb" + rb"nu0\d|/ho" + rb"me/jb" + rb"nu")
    hits = [f for f in tracked() if (REPO / f).is_file() and needle.search((REPO / f).read_bytes())]
    assert hits == []


def test_every_file_of_the_manifest_exists():
    rows = list(csv.DictReader(open(REPO / "MANIFEST.csv", encoding="utf-8")))
    assert rows and all((REPO / r["path"]).exists() for r in rows)


def test_ledger_numbers_agree_with_their_sources():
    r = subprocess.run([sys.executable, str(REPO / "reproduce" / "step1_check_numbers.py")], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-800:]


def test_the_english_documents_are_not_empty():
    docs = pathlib.Path(__file__).resolve().parents[1] / "docs"
    for name, least in (("fixed_analysis_document.md", 100_000), ("post_verdict_plan.md", 100_000), ("judge_specification.md", 20_000), ("verdict_record.md", 20_000), ("code_guide.md", 5_000), ("glossary.md", 1_000)):
        assert (docs / name).stat().st_size >= least, name


def test_no_text_file_is_empty_except_package_markers():
    root = pathlib.Path(__file__).resolve().parents[1]
    tracked = subprocess.run(["git", "ls-files"], cwd=root, capture_output=True, text=True, check=True).stdout.split("\n")[:-1]
    empty = [f for f in tracked if (root / f).stat().st_size == 0 and not f.endswith("__init__.py")]
    assert not empty, empty[:5]
