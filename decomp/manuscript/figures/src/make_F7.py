"""make_F7 draws Supplementary Fig. S2 (Figure 7 of earlier manuscript versions): the frozen judge on traveling heads, the acquisition matrix alone."""
import sys

sys.dont_write_bytecode = True            # no __pycache__ beside the scripts, in the ledger folder or in the repository

import datetime  # noqa: E402
import hashlib  # noqa: E402
import importlib.util  # noqa: E402
import io  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import pathlib  # noqa: E402
import subprocess  # noqa: E402

os.environ["CUDA_VISIBLE_DEVICES"] = ""   # CPU only

_HERE = pathlib.Path(__file__).resolve().parent
# v13: the judge-matrix script of v12 (v12 make_F6.py, copied as judge_matrix_v12.py) holds the shared reading, guards,
# caption check and save of the after-the-verdict figures; loaded by path (nothing is drawn or written on import).
_spec = importlib.util.spec_from_file_location("judge_matrix_shared", _HERE / "judge_matrix_v12.py")
B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(B)

fs, FC = B.fs, B.FC
REPO, WORK = B.REPO, B.WORK
DRAFTS = WORK / "stage1_work" / "drafts"
X14 = DRAFTS / "x14"
LEDGER_DIR = B.LEDGER_DIR
REMEASURE = WORK / "remeasure"
B2K = REMEASURE / "explore" / "b2000"
RULES_REC = DRAFTS / "x15" / "v12" / "RULES_RECORD.jsonl"
REREAD_LOG = X14 / "REREAD_LOG.jsonl"
OUT = pathlib.Path(B.TARGET)
SCRIPT = pathlib.Path(__file__).resolve()

RUN_ID = "F7-" + datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
WHAT = "v13d 그림 7 — 판단기 행렬만(Prisma 값은 표 2 · 장부) · 새 계산 없음 · 판정 불변 #34"
NAME = "F7_judge"
LEDGER_SRC_IDS = ("N-413", "N-116", "N-382", "N-416", "N-428", "N-029", "N-028", "N-423", "N-452", "N-493", "N-459")
READS7 = []                               # every file re-read for the figure (logged)


def stop(msg):
    sys.exit(f"STOP F7 (v13): {msg}")


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


def rel_repo(p):
    return os.path.relpath(pathlib.Path(p).resolve(), REPO)


def log_entry(log, pred):
    """(1-based line, entry) of the last entry of a JSONL log that satisfies pred."""
    lines = log.read_bytes().decode("utf-8").split("\n")
    for i in range(len(lines) - 1, -1, -1):
        s = lines[i].strip()
        if not s:
            continue
        try:
            e = json.loads(s)
        except json.JSONDecodeError:
            continue
        if pred(e):
            return i + 1, e
    stop(f"{log.name}: the entry just written was not found")


def reread(path, what, want_sha=None):
    """A derived file: logged with the ledger's reread_log.py (it hashes and appends under its lock), then read once; the
    bytes used are the bytes logged."""
    path = pathlib.Path(path).resolve()
    w = f"{WHAT} · {what} · 실행 {RUN_ID}"
    r = subprocess.run([sys.executable, "-I", "-B", str(LEDGER_DIR / "reread_log.py"), str(path), w], capture_output=True,
                       text=True)
    if r.returncode != 0:
        stop(f"reread_log.py {path}: {r.stdout.strip()} {r.stderr.strip()}")
    h = r.stdout.strip().splitlines()[-1].strip()
    raw = path.read_bytes()
    if sha_bytes(raw) != h:
        stop(f"{path} changed between the logged hash and the read")
    if want_sha and h != want_sha:
        stop(f"{path}: sha256 {h} is not the recorded {want_sha}")
    ln, e = log_entry(REREAD_LOG, lambda x: x.get("다시읽음") == str(path) and x.get("sha256") == h and x.get("무엇") == w)
    READS7.append(dict(kind="derived file, re-read logged by reread_log.py", file=str(path), sha256=h, bytes=len(raw),
                       log=f"drafts/x14/{REREAD_LOG.name}", line=ln, time=e.get("시각"), what=what, run=RUN_ID))
    return raw


def log_originals():
    """The repository files read_values() opens itself (judgements.json, cards.json) and the remeasure originals behind the
    ledger copies it reads, logged as re-read before the run; each must still equal its ledger copy."""
    seen = {}
    for i in LEDGER_SRC_IDS:
        r = B.LEDGER[i]
        rel, sha = r["근거 파일"], r["근거 파일 sha256"]
        if rel in seen:
            continue
        p = (REPO / rel).resolve()
        if not (str(p).startswith(str(REMEASURE.resolve()) + os.sep) or rel.startswith("decomp/results/")):
            continue
        reread(p, "장부 사본의 원본 — 판정 뒤 그림 값(7a · 7b · 7c 표지) · 장부 sha256 대조", want_sha=sha)
        seen[rel] = sha
    return seen


def draw7(v):
    """The figure is draw6()'s own (the former panel (a), without the panel letter)."""
    fig = B.draw6(v)
    B.parallel_ok(fig, v["_labels"], B.CS)
    return fig


def relabel(s):
    """v12 labels -> v13d Fig. 7: matrix «6» -> «7»; values guarded but not drawn here: Prisma «7a», «7b» -> Table 2, noise floor «7c» -> Fig. S5."""
    for old, new in (("6", "7"), ("7a", "guard (Table 2)"), ("7b", "guard (Table 2)"), ("7c", "guard (Fig. S5)")):
        if s.startswith(old):
            return new + s[len(old):]
    return s


def relabel_all():
    for e in B.DRAWN:
        e["panel"] = relabel(e["panel"])
    for e in B.SHOWN:
        e["where"] = relabel(e["where"])
    for i, roles in B.USED.items():
        B.USED[i] = [relabel(r) for r in roles]


def estimator_words(fig):
    """v17j (2026-11-20): the printed word «judge» becomes «estimator», as in the article, and the stale «Section 4.4» pointer
    becomes «S-M5»; nothing else changes.
    The overlap check of B.finish() runs afterwards, on the printed texts."""
    from matplotlib.text import Text
    n = 0
    for t in fig.findobj(Text):
        s = t.get_text()
        if "judge" in s:
            t.set_text(s.replace("the judge's own", "the estimator's own").replace("to the judge", "to the estimator"))
            n += 1
        if "Section 4.4; Supplementary Table S5" in t.get_text():          # stale since v13: the main text has no Section 4.4
            t.set_text(t.get_text().replace("Section 4.4; Supplementary Table S5", "S-M5; Table S5"))
    left = [t.get_text() for t in fig.findobj(Text) if "judge" in t.get_text().lower()]
    if n == 0 or left:
        raise SystemExit(f"STOP F7: estimator words: replaced {n}, left {left[:3]}")


def reads_md(head):
    L = [head, "", f"Run `{RUN_ID}` ({now()}). No image, mask, transform or template was read by this script.", "",
         "| # | kind | file | sha256 | bytes | log : line | time | what |", "|---|---|---|---|---|---|---|---|"]
    for i, x in enumerate(READS7, 1):
        f = x["file"] if not str(x["file"]).startswith("/") else rel_repo(x["file"])
        L.append(f"| {i} | {x['kind']} | `{f}` | `{x['sha256']}` | {x['bytes']} | {x['log']}:{x['line']} | {x.get('time')} | "
                 f"{x['what']} |")
    L += ["", "Ledger copies read by judge_matrix_v12.read_values() (drafts/x15/ledger/sources/<file>/<sha12>; checked "
          "against the ledger sha256, not logged — they are the ledger's own copies):", ""]
    for rel, sha in sorted(B.READS.items()):
        L.append(f"- `{rel}` `{sha[:12]}`")
    return "\n".join(L) + "\n"


def merged_srcs(*dicts):
    out, seen = {}, set()
    for d in dicts:
        for k, s in d.items():
            if s.rel not in seen:
                seen.add(s.rel)
                out[f"{len(out)}_{k}"] = s
    return out


SRC7 = ("ledger copies of `decomp/results/remeasure/judgements.json` `V4.cells` · `V4.summary`, "
        "`remeasure/explore/d1/d1_rulethr.json` `cells`, `remeasure/explore/d1/d1_post_sym.json` `sym`, "
        "`remeasure/explore/rv/rv_4a_coil.json` `by_device` · `sessions`, `decomp/results/stage3/cards/cards.json`, "
        "`decomp/s4/plan.py:17–22` `ONH`, `decomp/stage1/analyze3.py` `NOT2ING_FIXED` (drawn); "
        "`remeasure/explore/cause/cause1b_allunits.json` `sets.all`, `remeasure/explore/cause/cause1_raw.json` `expected`, "
        "`remeasure/explore/suff/suff_c3.json` `cells`, `decomp/s4/s4cfg.py` `COIL_L` (guarded, not drawn) · drawn by "
        "`decomp/manuscript/figures/src/make_F7.py` with the reading, guards and matrix drawing of `judge_matrix_v12.py`, "
        "which stop unless every printed number equals its ledger row and every drawn value equals its field in the "
        "ledger copies and, where one exists, its ledger row")
# Caption: direction S50 A5.
CAPTION = ("![**Fig. 7. The frozen estimator on traveling heads.** Acquisitions ordered by TE, with vendor, TE, TR, receive coil and "
           "direction scheme. Cells: coverage of zero minus the cell's fixed threshold, {N-380} ODI cells; upper triangle, real data "
           "(blinded); lower, the estimator's own simulated test data scored with the real-data rule (after the verdict; Supplementary "
           "Table S5). Check mark, threshold met; outlined, the Prisma pair (Table 2); below the matrix, its NDI cell (unblinded "
           "check); hatched, no data.](figures/F7_judge.pdf){#fig:fig7 .wide width=100%}")


def notes7(v):
    return [
        "v13d (direction S50 A5, author decision 2026-11-08): the figure is v12 Fig. 6 (the acquisition matrix), drawn by "
        "judge_matrix_v12.draw6() unchanged. The former panels (b) and (c) became main-text Table 2 (ledger values); v12 Fig. 7c "
        "moved to Supplementary Fig. S5 in v13.",
        "Sources and guards (formerly the [src] blocks of the figure lines): " + SRC7 + ".",
    ] + B.notes6(v)


def main():
    log_originals()
    v = B.read_values()                                   # every data guard of judge_matrix_v12.py; nothing written before
    fig = draw7(v)
    relabel_all()
    estimator_words(fig)
    extra = ["## Layout", "", f"Figure 190 × {round(fig.get_size_inches()[1] * 25.4)} mm: judge_matrix_v12.draw6(), its own height.", ""]
    srcs = merged_srcs(v["srcs_judge"], v["srcs7"])
    audit = B.finish(fig, 7, NAME, CAPTION, srcs, notes7(v), extra)
    (OUT / "F7_reads.md").write_text(reads_md("# v13 Fig. 7 — files read for the figure (by make_F7.py)"), encoding="utf-8")
    (OUT / "F7_reads.json").write_text(json.dumps(READS7, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    cvd = FC.cvd_sheet(OUT / f"{NAME}.png", OUT / "cvd" / "F7_cvd.png")
    print("cvd sheet:", cvd.name, hashlib.sha256(cvd.read_bytes()).hexdigest()[:16])
    print("files re-read (logged):", len(READS7), "· ledger copies read:", len(B.READS))
    print(audit)


if __name__ == "__main__":
    main()
