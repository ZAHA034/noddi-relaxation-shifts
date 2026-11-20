"""Draws Figure 6 of manuscript v12: the frozen judge on traveling heads as one acquisition matrix, from stored values only."""
import ast
import csv
import hashlib
import importlib.util
import json
import math
import os
import pathlib
import re
import subprocess
import sys
from decimal import Decimal, ROUND_HALF_UP

sys.dont_write_bytecode = True        # never leave __pycache__ beside the scripts or in the ledger folder
_FILE = pathlib.Path(__file__).resolve()
_HERE = _FILE.parent


def _find_repo():
    """DMRI_REPO, else the repository this file sits in (decomp/manuscript/figures/src), else the one beside data."""
    env = os.environ.get("DMRI_REPO")
    if env:
        return pathlib.Path(env)
    if len(_FILE.parents) > 4 and (_FILE.parents[4] / "decomp" / "manuscript").is_dir():
        return _FILE.parents[4]
    return _FILE.parents[6].parent / "dmri_cov"          # <project-notes>/x15/ms/S41_files/src


REPO = _find_repo()
WORK = pathlib.Path(os.environ.get("DMRI_WORK") or REPO / "data")
if not (REPO / "decomp" / "manuscript" / "figures" / "src").is_dir() or not (WORK / "stage1_work").is_dir():
    sys.exit(f"STOP F6/F7: repository {REPO} or work folder {WORK} not found (set DMRI_REPO / DMRI_WORK)")
if not os.environ.get("FIG_OUT"):
    sys.exit("STOP F6/F7: set FIG_OUT (the only folder this script writes to)")
STYLE = str(REPO / "decomp" / "manuscript" / "figures" / "src")
sys.path.insert(0, STYLE)
import figstyle as fs  # noqa: E402

TARGET = str(pathlib.Path(os.environ["FIG_OUT"]).resolve())
pathlib.Path(TARGET).mkdir(parents=True, exist_ok=True)
fs.OUTDIR = TARGET
import fig_data_templates_v3 as T  # noqa: E402  acquisition types, no-data cells, pair names, vendors; nothing drawn
import matplotlib  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402
import matplotlib.text as mtext  # noqa: E402

_spec = importlib.util.spec_from_file_location("figcolors", _HERE / "figcolors.py")
FC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(FC)          # S42 §1.2 colours; asserts the palette rules on import

LEDGER_DIR = WORK / "stage1_work" / "drafts" / "x15" / "ledger"
sys.path.insert(0, str(LEDGER_DIR))
import ledger_check as LC  # noqa: E402  (eval_field · same_value · rounding_ok — the ledger's own semantics)

OUT = pathlib.Path(TARGET)
NAME = "F6_judge_matrix"
SRC = LEDGER_DIR / "sources"
LEDGER = {r["ID"]: r for r in csv.DictReader((LEDGER_DIR / "numbers.csv").open(encoding="utf-8"))}
SUPP_MD = pathlib.Path(os.environ.get("SUPP_MD") or WORK / "stage1_work" / "drafts" / "x15" / "ms" / "S40_files" /
                       "supplement.md")
S5_MODE = os.environ.get("S5_MODE", "part_e")

INK, INK2, MUTED, GRID, AXIS, SURF = fs.INK, fs.INK2, fs.MUTED, fs.GRID, fs.AXIS, fs.SURF
PT, PTL = fs.PT_TEXT, fs.PT_LETTER
ZEBRA = "#f6f6f3"
SIM = fs.GREYS[1]                      # simulator / after-the-verdict marks (contrast 6.6 : 1 on white) — used by Fig. 7
SIM_LIGHT = fs.GREYS[2]                # recalculated thresholds, Philips/GE (3.6 : 1) — used by Fig. 7
MRK = {"NDI": "o", "ODI": "^"}         # coordinates by shape (black), not by colour
MS = {"NDI": 3.9, "ODI": 4.3}
GIVEN = {"OXF1PRI": "32", "OXF2PRI": "64"}       # receive coils given to the judge (asserted from plan.py and the candidates)
CAND_ALL = ("32", "64")
TYPE = {"match": 0, "te": 1, "tr": 2, "tetr": 3, "tetr_design": 4}     # judge_v4 acq_type -> old Fig. 4 type (asserted)
SIEMENS = ("OXF1PRI", "OXF2PRI", "OXF3TRI")
OTHERS = ("NOT1ACH", "NOT2ING-89", "NOT2ING-98", "NOT4GEP", "OXF4GEP")
NO_COIL = ("NOT1ACH", "NOT2ING", "NOT4GEP", "OXF4GEP")                # sidecars without a receive-coil record (asserted)
MODELS = ("primary", "s1_inv_snr", "s2_window_intercept")             # the three regression models (s3 = primary, asserted)
SEP = "~"                                                             # pair names as in the text
TAG_SOLID, TAG_DASHED = "-", (0, (2.2, 1.6))                          # fixed tests: solid border; after the verdict: dashed
SCHEME = {"S": "A", "V": "B"}                                         # plan.py ONH design letters -> scheme letters shown


def stop(msg):
    sys.exit(f"STOP F6/F7: {msg}")


def pname(p):
    return T.pname(p, sep=SEP)


# ----------------------------------------------------------------------------- ledger
KEYRX = re.compile(r'\["([^"]*)"\]|\[\'([^\']*)\'\]|\[(\d+)\]|([^.\[\]]+)')


def keys_of(path):
    out = []
    for m in KEYRX.finditer(path):
        q1, q2, idx, key = m.groups()
        out.append(int(idx) if idx is not None else (q1 if q1 is not None else q2 if q2 is not None else key))
    return tuple(out)


FIELD_INDEX = {}                       # (file, sha256, keys) -> ledger IDs whose 근거 필드 (basis field) is that json field
for _i, _r in LEDGER.items():
    _kind, _, _rest = _r["근거 필드"].partition(":")
    if _kind == "json":
        FIELD_INDEX.setdefault((_r["근거 파일"], _r["근거 파일 sha256"], keys_of(_rest)), []).append(_i)

ROWCHECK = {}                          # ID -> re-derived value (None for lit · none)
USED = {}                              # ID -> list of roles («6…» or «7…» first, for the figure's values file)
DRAWN = []                             # every drawn value: panel, element, file, sha, field, value, ledger IDs
SHOWN = []                             # every number printed in a figure: where, text, ledger IDs
TABLED = []                            # every value this script prints in Supplementary Table S5 (after the verdict)
READS = {}                             # ledger copies read: file -> sha256 (Src objects and re-derived ledger rows)


def check(id_, value=None, role=""):
    """Re-derive ledger row id_ from its source copy, re-check its rounding, compare a drawn value; return its 원고 표기 (manuscript notation)."""
    if id_ not in LEDGER:
        stop(f"{id_} is not in the ledger")
    r = LEDGER[id_]
    if id_ not in ROWCHECK:
        if r["근거 파일"] and not r["근거 필드"].startswith(("lit:", "none:")):
            p = SRC / r["근거 파일"].replace("../", "__up__/") / r["근거 파일 sha256"][:12]
            if hashlib.sha256(p.read_bytes()).hexdigest() != r["근거 파일 sha256"]:
                stop(f"{id_}: ledger copy {p} differs from the ledger sha256")
            READS[r["근거 파일"]] = r["근거 파일 sha256"]
        v = LC.eval_field(SRC, r)
        if v is not None and not LC.same_value(v, r["원 값"]):
            stop(f"{id_}: source gives {v!r}, ledger raw value {r['원 값']!r}")
        ok, why = LC.rounding_ok(r["원 값"], r["원고 표기"], r["반올림 규칙"])
        if ok is False:
            stop(f"{id_}: 원고 표기 {r['원고 표기']!r} is not the raw value rounded by «{r['반올림 규칙']}» ({why})")
        ROWCHECK[id_] = v
    if value is not None and not LC.same_value(value, r["원 값"]):
        stop(f"{id_}: drawn value {value!r} differs from the ledger raw value {r['원 값']!r}")
    if role:
        USED.setdefault(id_, []).append(role)
    return r["원고 표기"]


def show(where, ids, text):
    """Register a printed number: text built from 원고 표기 (manuscript notation) only."""
    if not ids:
        stop(f"{where}: a printed number needs a ledger row")
    SHOWN.append(dict(where=where, text=text, ids=list(ids)))
    return text


class Src:
    """The ledger's frozen copy of one source file, chosen through a ledger row that cites it."""

    def __init__(self, via_id):
        r = LEDGER[via_id]
        self.rel, self.sha = r["근거 파일"], r["근거 파일 sha256"]
        self.path = SRC / self.rel.replace("../", "__up__/") / self.sha[:12]
        self.raw = self.path.read_bytes()
        if hashlib.sha256(self.raw).hexdigest() != self.sha:
            stop(f"{self.path} differs from the ledger sha256 {self.sha}")
        READS[self.rel] = self.sha
        self.doc = json.loads(self.raw) if self.rel.endswith(".json") else None
        self.text = self.raw.decode("utf-8")

    def get(self, *keys):
        cur = self.doc
        for k in keys:
            cur = cur[k]
        return cur


def fmt_keys(keys):
    out = ""
    for k in keys:
        out += f"[{k}]" if isinstance(k, int) else (f'["{k}"]' if "." in str(k) else ("." if out else "") + str(k))
    return out


def drawn(panel, element, src, *keys):
    """Read one stored value, compare it with every ledger row that cites the same file (sha256) and field, log it."""
    v = src.get(*keys)
    ids = FIELD_INDEX.get((src.rel, src.sha, tuple(keys)), [])
    for i in ids:
        check(i, v, f"{panel}: {element}")
    DRAWN.append(dict(panel=panel, element=element, file=src.rel, sha=src.sha[:12], field="json:" + fmt_keys(keys),
                      value=v, ids=list(ids)))
    return v


def drawn_md(panel, element, src, line, value, ids):
    DRAWN.append(dict(panel=panel, element=element, file=src.rel, sha=src.sha[:12], field=f"md:{line}", value=value,
                      ids=list(ids)))


def drawn_py(panel, element, src, field, value, ids=()):
    """A drawn value taken from a literal assignment of a Python source copy (read with ast; nothing executed)."""
    DRAWN.append(dict(panel=panel, element=element, file=src.rel, sha=src.sha[:12], field="py:" + field, value=value,
                      ids=list(ids)))
    return value


def tabled(element, src, keys, value, printed):
    """A value this script prints in Supplementary Table S5: compared with any ledger row on the same field, logged."""
    ids = FIELD_INDEX.get((src.rel, src.sha, tuple(keys)), [])
    for i in ids:
        check(i, value, f"6 (Table S5): {element}")
    TABLED.append(dict(element=element, file=src.rel, sha=src.sha[:12], field="json:" + fmt_keys(keys), value=value,
                       printed=printed, ids=list(ids)))


def same(a, b, what):
    if a is None or b is None or abs(float(a) - float(b)) > 1e-12:
        stop(f"{what}: {a} != {b}")


def py_assign(src, name):
    """A top-level literal assignment in a Python source copy, read with ast.literal_eval: (value, first, last line)."""
    tree = ast.parse(src.text)
    hits = [n for n in tree.body if isinstance(n, ast.Assign) and len(n.targets) == 1
            and isinstance(n.targets[0], ast.Name) and n.targets[0].id == name]
    if len(hits) != 1:
        stop(f"{src.rel}: {len(hits)} top-level assignments of {name}")
    return ast.literal_eval(hits[0].value), hits[0].lineno, hits[0].end_lineno


def half_up(x, k):
    """The ledger's rounding (ledger_check._half_up): half up at k decimals, from 12 significant digits."""
    return Decimal(f"{x:.12g}").quantize(Decimal(1).scaleb(-k), rounding=ROUND_HALF_UP)


S5_TIES = []                           # Table S5 entries printed half to even at an exact binary tie (ledger rule: half up)
S5_NEW_TIES = []                       # the after-the-verdict values this script prints at an exact tie (printed half up)


def printed_is(text, value, signed=False, where=""):
    """«text» is «value» rounded to the decimals printed (with «+» / «−» when signed and not zero): «half up» (the
    ledger's rule), «half to even» (Python's format at an exact binary tie; recorded in S5_TIES), or False."""
    t = text.strip().replace("−", "-")
    if not re.fullmatch(r"[+-]?\d+(?:\.\d+)?", t):
        return False
    k = len(t.split(".")[1]) if "." in t else 0
    want, py = half_up(value, k), Decimal(f"{value:.{k}f}")
    if Decimal(t) == want:
        how = "half up"
    elif Decimal(t) == py:
        how = "half to even"
        S5_TIES.append(f"{where}: {text} = {value!r} rounded half to even (half up would print {want})")
    else:
        return False
    if signed:
        ok = t.startswith("+") if Decimal(t) > 0 else t.startswith("-") if Decimal(t) < 0 else t[0].isdigit()
    else:
        ok = not t.startswith(("+", "-")) or Decimal(t) < 0
    return how if ok else False


# ----------------------------------------------------------------------------- Table S5 of the supplement (read only)
S5_HEAD = ["Cell", "Card pair", "ΔTE (ms)", "ΔTR (s)", "Same scheme", "Coils given", "Comp.", "Role", "Subjects",
           "Units: before → used", "Excluded: misspec. / k unfilled", "Coverage", "Threshold", "Result", "Gap flag",
           "Fixed sentence"]
S5B_HEAD = ["Cell", "Card pair · component", "By coil combination (A, B): units over the threshold / units lacking samples"]


def md_cells(line):
    s = line.strip()
    if not (s.startswith("|") and s.endswith("|")):
        return None
    return [c.strip() for c in re.split(r"(?<!\\)\|", s[1:-1])]


def parse_supp(path):
    raw = path.read_bytes()
    text = raw.decode("utf-8")
    L = text.split("\n")

    def find(pred, start, what):
        for i in range(start, len(L)):
            if pred(L[i]):
                return i
        stop(f"supplement {path.name}: {what} not found")

    h = find(lambda s: s.strip() == "### Table S5 {.unnumbered}", 0, "«### Table S5 {.unnumbered}»")
    cap = find(lambda s: s.startswith("Table: **Table S5.**"), h, "the Table S5 caption")
    head = find(lambda s: s.startswith("| Cell | Card pair |"), cap, "the Table S5 header")
    if md_cells(L[head]) != S5_HEAD:
        stop(f"Table S5 header {md_cells(L[head])}")
    if not re.fullmatch(r"\|(?:-+\|)+", L[head + 1].strip()) or md_cells(L[head + 1]) is None \
            or len(md_cells(L[head + 1])) != len(S5_HEAD):
        stop("Table S5 separator line")
    rows, i = [], head + 2
    while i < len(L) and L[i].startswith("|"):
        c = md_cells(L[i])
        if c is None or len(c) != len(S5_HEAD):
            stop(f"Table S5 line {i + 1}: {len(c or [])} cells")
        rows.append(dict(zip(S5_HEAD, c), _line=i + 1, _text=L[i]))
        i += 1
    last = i - 1
    notes = find(lambda s: s.startswith("*Notes.*"), i, "the Table S5 notes")
    bcap = find(lambda s: s.startswith("Table: **Table S5 (b).**"), i, "the Table S5 (b) caption")
    bhead = find(lambda s: s.startswith("| Cell | Card pair · component |"), bcap, "the Table S5 (b) header")
    if md_cells(L[bhead]) != S5B_HEAD:
        stop(f"Table S5 (b) header {md_cells(L[bhead])}")
    brows, j = {}, bhead + 2
    while j < len(L) and L[j].startswith("|"):
        c = md_cells(L[j])
        brows[c[0]] = dict(pair_comp=c[1], combos=re.findall(r"\((\d+), (\d+)\):", c[2]), _line=j + 1)
        j += 1
    s6 = [k for k, s in enumerate(L) if s.strip() == "### Table S6 {.unnumbered}"]
    if len(s6) != 1:
        stop("«### Table S6 {.unnumbered}» is not unique")
    return dict(path=path, sha=hashlib.sha256(raw).hexdigest(), text=text, lines=L, head=head, last=last, cap=cap,
                notes=notes, rows=rows, b=brows, s6=s6[0])


# ----------------------------------------------------------------------------- values (Fig. 6 and Fig. 7)
def read_values():
    J = Src("N-413")                    # decomp/results/remeasure/judgements.json
    repo_j = REPO / "decomp" / "results" / "remeasure" / "judgements.json"
    if hashlib.sha256(repo_j.read_bytes()).hexdigest() != J.sha:
        stop("the repository judgements.json differs from the ledger copy")
    D1 = Src("N-382")                   # remeasure/explore/d1/d1_rulethr.json
    PS = Src("N-416")                   # remeasure/explore/d1/d1_post_sym.json
    RV = Src("N-428")                   # remeasure/explore/rv/rv_4a_coil.json (sidecar coil fields, read after the verdict)
    CD = Src("N-116")                   # decomp/results/stage3/cards/cards.json (dTE · dTR_ms · design_differs)
    PLN = Src("N-250")                  # decomp/s4/plan.py (ONH: scheme letter · TE range ms · TR range s · coil given)
    A3 = Src("N-114")                   # decomp/stage1/analyze3.py (NOT2ING_FIXED: the confirmed NOT2ING TR ranges, s)
    C1 = Src("N-029")                   # remeasure/explore/cause/cause1_raw.json
    C1B = Src("N-028")                  # remeasure/explore/cause/cause1b_allunits.json
    SU = Src("N-423")                   # remeasure/explore/suff/suff_c3.json
    BS = Src("N-452")                   # remeasure/explore/b2000/b2k_result_siemens.json
    BA = Src("N-493")                   # remeasure/explore/b2000/b2k_result_all.json
    KS = Src("N-459")                   # remeasure/explore/b2000/kat_summary_siemens.json
    PL = Src("N-462")                   # drafts/EXPLORATORY_PLAN_2026-11-01.md (copy 5059d6d1cb55)
    CF = Src("N-429")                   # decomp/s4/s4cfg.py
    repo_cd = REPO / "decomp" / "results" / "stage3" / "cards" / "cards.json"
    if hashlib.sha256(repo_cd.read_bytes()).hexdigest() != CD.sha:
        stop("the repository cards.json differs from the ledger copy")
    srcs_judge = dict(J=J, D1=D1, PS=PS, RV=RV, CD=CD, PLN=PLN, A3=A3)
    srcs7 = dict(J=J, PS=PS, C1=C1, C1B=C1B, SU=SU, BS=BS, BA=BA, KS=KS, PL=PL, CF=CF)
    V4 = J.get("V4")
    cj, summ = V4["cells"], V4["summary"]
    d1c = D1.get("cells")

    # ---------------- Fig. 6 tracks: what the sidecars record and what the judge was given
    byd = RV.get("by_device")
    for dev, ch in (("OXF1PRI", "32"), ("OXF2PRI", "64"), ("OXF3TRI", "32")):
        if byd[dev]["channels"] != {ch: byd[dev]["n_sessions"]}:
            stop(f"{dev}: sidecar channels {byd[dev]['channels']} are not {ch} in every session")
    if byd["OXF3TRI"]["ReceiveCoilName"] != {"32Ch_Head": byd["OXF3TRI"]["n_sessions"]}:
        stop("OXF3TRI: the sidecars do not all record 32Ch_Head")
    for dev in NO_COIL:
        if byd[dev]["channels"] != {"None": byd[dev]["n_sessions"]}:
            stop(f"{dev}: a receive coil is recorded, so «not given» would need another reason")
    if set(byd) != set(SIEMENS) | set(NO_COIL):
        stop(f"rv_4a_coil devices {sorted(byd)}")
    check("N-428", RV.get("sessions", "sub-03286_ses-OXF3TRI001", "channels"),
          "6 (guard): OXF3TRI sidecar coil, 32 channels, not given to the judge (open chip)")

    # ---------------- Fig. 6 tracks: TE, TR, scheme and coil of the eight acquisitions (stage-3 card values)
    onh, l0, l1 = py_assign(PLN, "ONH")
    if (l0, l1) != (17, 22):
        stop(f"plan.py ONH is at lines {l0}–{l1}, not 17–22 as Table S5 cites")
    nfix, f0, f1 = py_assign(A3, "NOT2ING_FIXED")
    if set(nfix) != {"NOT2ING-89", "NOT2ING-98"}:
        stop(f"NOT2ING_FIXED keys {sorted(nfix)}")
    if set(onh) != set(SIEMENS) | set(OTHERS):
        stop(f"plan.py ONH acquisitions {sorted(onh)}")
    if {x[0] for x in onh.values()} != set(SCHEME):
        stop(f"plan.py ONH scheme letters {sorted({x[0] for x in onh.values()})}")
    acq = {}
    for nm, (letter, te_r, tr_r, coil) in onh.items():
        tr_src = nfix.get(nm, tr_r)
        te, tr = (te_r[0] + te_r[1]) / 2, (tr_src[0] + tr_src[1]) / 2           # as analyze3.py x11 (cards)
        dev = nm.split("-")[0]
        ch = byd[dev]["channels"]
        rec = None if list(ch) == ["None"] else int(list(ch)[0])
        given = None if coil is None else int(coil)
        if (str(given) if given else None) != GIVEN.get(nm):
            stop(f"{nm}: plan.py gives the judge coil {coil}, not {GIVEN.get(nm)}")
        if given is not None and rec != given:
            stop(f"{nm}: coil given {given} but sidecars record {rec}")
        acq[nm] = dict(vendor=T.vendor(nm), te=te, tr=tr, letter=letter, scheme=SCHEME[letter], given=given, rec=rec,
                       chip="filled" if given is not None else "open" if rec is not None else "none")
        el = f"track {nm}"
        drawn_py("6 tracks", f"{el}: TE (dot), midpoint of the recorded range", PLN, f"ONH['{nm}'][1] = {list(te_r)}", te)
        if nm in nfix:
            drawn_py("6 tracks", f"{el}: TR (dot, s), midpoint of the confirmed range (as the cards)", A3,
                     f"NOT2ING_FIXED['{nm}'] = {list(tr_src)}", tr)
        else:
            drawn_py("6 tracks", f"{el}: TR (dot, s), midpoint of the recorded range", PLN, f"ONH['{nm}'][2] = {list(tr_r)}",
                     tr)
        drawn_py("6 tracks", f"{el}: scheme letter {SCHEME[letter]} (from «{letter}»)", PLN, f"ONH['{nm}'][0]", letter)
        drawn_py("6 tracks", f"{el}: receive coil given to the judge (chip filled when given)", PLN, f"ONH['{nm}'][3]",
                 given)
        DRAWN.append(dict(panel="6 tracks", element=f"{el}: receive coil recorded in the sidecars (chip open when not given)",
                          file=RV.rel, sha=RV.sha[:12], field=f"json:by_device.{dev}.channels", value=ch, ids=[]))
    order = sorted(acq, key=lambda k: (acq[k]["te"], k))
    cards = [c for c in CD.doc if c.get("set") == "ON-Harmony"]
    if len(cards) != 28 or sorted(tuple(c["pair"].split("~")) for c in cards) != sorted(T.DESIGN):
        stop("cards.json ON-Harmony pairs are not the 28 card pairs")
    for c in cards:
        a, b = c["pair"].split("~")
        if abs((acq[b]["te"] - acq[a]["te"]) - c["dTE"]) > 1e-9:
            stop(f"{c['pair']}: track TE difference {acq[b]['te'] - acq[a]['te']} != cards dTE {c['dTE']}")
        if abs(1000 * (acq[b]["tr"] - acq[a]["tr"]) - c["dTR_ms"]) > 1e-6:
            stop(f"{c['pair']}: track TR difference != cards dTR_ms {c['dTR_ms']}")
        if (acq[a]["letter"] != acq[b]["letter"]) != c["design_differs"]:
            stop(f"{c['pair']}: scheme letters disagree with cards design_differs")
    check("N-116", sum(c["design_differs"] for c in cards), "6 (guard): card pairs with different schemes (scheme track)")
    tes, trs = [acq[k]["te"] for k in acq], [acq[k]["tr"] for k in acq]
    check("N-250", [min(tes), max(tes)], "6 (guard): TE range of the drawn TE track (68–98 ms)")
    check("N-630", 1000 * min(trs), "6 (guard): shortest TR of the drawn TR track (3,600 ms)")
    check("N-631", 1000 * max(trs), "6 (guard): longest TR of the drawn TR track (10,215 ms)")

    # ---------------- Fig. 6: cells
    with_data = [k for k, c in cj.items() if c["verdict"] != "자료 없음"]
    nodata_keys = [k for k, c in cj.items() if c["verdict"] == "자료 없음"]
    if sorted(with_data) != sorted(d1c):
        stop("d1_rulethr cells differ from the calibration-test cells with data")
    nd_pairs = [tuple(k.split("|")[0].split("~")) for k in nodata_keys]
    if sorted(nd_pairs) != sorted(T.NO_DATA) or not all(k.endswith("|ODI") for k in nodata_keys):
        stop(f"no-data cells {nd_pairs}")
    nodata = sorted(T.NO_DATA, key=T.acq_type)
    cells = []
    for k in with_data:
        c, d = cj[k], d1c[k]
        pair_s, comp = k.split("|")
        p = tuple(pair_s.split("~"))
        if p not in T.DESIGN or TYPE[c["acq_type"]] != T.acq_type(p) or d["acq_type"] != c["acq_type"]:
            stop(f"{k}: acquisition type {c['acq_type']} disagrees with old Fig. 4")
        if (comp, c["blind"]) not in (("NDI", "비맹검 점검"), ("ODI", "맹검 확인")):
            stop(f"{k}: blinding {c['blind']}")
        if c["threshold_fallback"] or d["threshold_fallback_stored"]:
            stop(f"{k}: a fallback threshold would need its own mark")
        # ledger question (i): one fixed threshold per cell (a number, not one per k), while coverage varies over k
        if not isinstance(c["threshold"], (int, float)) or isinstance(c["threshold"], bool):
            stop(f"{k}: the fixed threshold is not one number per cell")
        el = f"cell {pair_s} {comp}"
        up = "6 upper" if comp == "ODI" else "6 NDI footnote"
        lo = "6 lower" if comp == "ODI" else "6 NDI footnote"
        cov = drawn(up, f"{el}: real coverage, median over the nine k (colour of the real-data cell, as coverage − fixed "
                    "threshold)", J, "V4", "cells", k, "coverage")
        cmin = J.get("V4", "cells", k, "coverage_min")         # guard only (the matrix does not draw the range over k)
        cmax = J.get("V4", "cells", k, "coverage_max")
        thr = drawn("6 both" if comp == "ODI" else "6 NDI footnote", f"{el}: fixed threshold (subtracted in both cells; "
                    "✓ at or above)", J, "V4", "cells", k, "threshold")
        tct = d["test_coverage_true_combo_kmed"]              # guard only (not drawn)
        tcr = drawn(lo, f"{el}: test-data coverage under the real-data rule, median over k (colour of the after-the-verdict "
                    "cell, as coverage − fixed threshold)", D1, "cells", k, "test_coverage_rule_kmed")
        thr_r = D1.get("cells", k, "threshold_rule_ours")     # not drawn: printed in Supplementary Table S5 (below)
        same(d["coverage_real"], cov, f"{k} real coverage (d1_rulethr vs judgements)")
        same(d["threshold_stored"], thr, f"{k} fixed threshold (d1_rulethr vs judgements)")
        if not (cmin <= cov <= cmax):
            stop(f"{k}: median outside the k range")
        if bool(c["passed"]) != (cov >= thr) or c["verdict"] != ("통과" if c["passed"] else "실패"):
            stop(f"{k}: stored verdict {c['verdict']} at coverage {cov} and threshold {thr}")
        opts_a = [GIVEN[p[0]]] if p[0] in GIVEN else list(CAND_ALL)
        opts_b = [GIVEN[p[1]]] if p[1] in GIVEN else list(CAND_ALL)
        cand = sorted(f"{x}|{y}" for x in opts_a for y in opts_b)
        if sorted(d["candidates"]) != cand:
            stop(f"{k}: candidates {d['candidates']} do not follow from the coils given to the judge")
        n = len(cand)
        if n == 2 and not (p[0] in GIVEN and p[1] not in GIVEN):
            stop(f"{k}: two candidates but not «A given, B not»")
        if n == 1:
            same(tct, tcr, f"{k}: with one candidate the two scoring rules coincide")
        if d["n_clusters"] != c["n_subjects"]:
            stop(f"{k}: threshold drawn with {d['n_clusters']} clusters, not the cell's {c['n_subjects']} participants")
        if c["n_units_used"] != c["n_units"] - c["n_excluded_misspec"] - c["n_units_dropped_k"]:
            stop(f"{k}: units used {c['n_units_used']} are not the units less the misspecified and k-incomplete ones")
        cells.append(dict(key=k, pair=p, comp=comp, acq=c["acq_type"], cov=cov, cmin=cmin, cmax=cmax, thr=thr,
                          passed=bool(c["passed"]), tct=tct, tcr=tcr, thr_r=thr_r, n_cand=n, cands=tuple(cand),
                          n_subj=c["n_subjects"], blind=c["blind"]))
    cells.sort(key=lambda e: ({1: 0, 2: 1, 4: 2}[e["n_cand"]], 0 if e["comp"] == "NDI" else 1, TYPE[e["acq"]],
                              pname(e["pair"])))
    shared = {}
    for e in cells:
        shared.setdefault((e["acq"], e["cands"], e["comp"]), []).append(e)
    for key, es in shared.items():
        for fld in ("tct", "tcr"):
            if max(x[fld] for x in es) - min(x[fld] for x in es) > 1e-12:
                stop(f"{key}: {fld} is not shared by the cells of one acquisition type and candidate set")
    ndi = [e for e in cells if e["comp"] == "NDI"]
    odi = [e for e in cells if e["comp"] == "ODI"]
    by_n = {n: [e for e in cells if e["n_cand"] == n] for n in (1, 2, 4)}
    if [e["pair"] for e in by_n[1]] != [T.REF_PAIR, T.REF_PAIR] or [e["comp"] for e in by_n[1]] != ["NDI", "ODI"]:
        stop("the group with both coils given must be OXF1PRI~OXF2PRI, NDI and ODI")
    # ledger question (ii): within the coil groups of two and four candidates, the fixed and the recalculated thresholds
    # differ between cells, and the applied-rule coverage takes more than one value -> one value per cell
    ans_ii = {}
    for n in (2, 4):
        g = by_n[n]
        ans_ii[n] = dict(thr={e["thr"] for e in g}, thr_r={e["thr_r"] for e in g}, tcr={e["tcr"] for e in g})
        if len(ans_ii[n]["thr"]) < 2 or len(ans_ii[n]["tcr"]) < 2:
            stop(f"group of {n} candidates: values shared after all; revisit the design")
    if len(ans_ii[2]["thr_r"]) < 2:
        stop("recalculated thresholds of the two-candidate cells are shared; revisit the design")
    # the matrix: every pair of the eight acquisitions once, as an ODI cell with data or a cell without data
    odi_pairs = [frozenset(e["pair"]) for e in odi] + [frozenset(p) for p in nodata]
    allp = {frozenset((a, b)) for i, a in enumerate(order) for b in order[i + 1:]}
    if len(odi_pairs) != len(set(odi_pairs)) or set(odi_pairs) != allp or len(allp) != 28:
        stop("the ODI cells with and without data are not the 28 pairs of the eight acquisitions, once each")
    # counts (guards; printed: 0 of 23, 22 of 24; caption: 23)
    check("N-021", len(cells), "6 lower tag line: cells with data (printed «22 of 24»)")
    check("N-149", len(cells), "6 (guard): cells with data, NDI + ODI")
    check("N-381", len(ndi), "6 (guard): NDI cells with data")
    check("N-380", len(odi), "6 (caption): ODI cells with data")
    check("N-255", len(nodata), "6 (guard): ODI cells without data (hatched)")
    check("N-266", sum(e["passed"] for e in odi), "6 upper tag line: ODI cells passed (printed «0 of 23»)")
    check("N-379", sum(e["passed"] for e in odi), "6 (guard): ODI cells passed")
    check("N-267", sum(not e["passed"] for e in odi), "6 (guard): ODI cells failed")
    check("N-268", sum(not e["passed"] for e in ndi), "6 (guard): NDI cell failed (printed as «NDI failed»)")
    check("N-020", sum(not e["passed"] for e in cells), "6 (guard): cells failed")
    check("N-378", sum(not e["passed"] for e in cells), "6 (guard): cells failed")
    for comp, n_ in (("NDI", len(ndi)), ("ODI", len(odi))):
        if (summ[comp]["with_data"], summ[comp]["passed"]) != (n_, sum(e["passed"] for e in cells if e["comp"] == comp)):
            stop(f"calibration-test summary {comp}")
    if summ["ODI"]["no_data"] != len(nodata) or summ["FISO"]["with_data"] != 0:
        stop("calibration-test summary no-data / FISO")
    # ranges given in the text of 4.6 (guards; the per-cell values are printed in Supplementary Table S5)
    g2, g4 = by_n[2], by_n[4]
    check("N-382", min(e["tcr"] for e in g2), "6 (guard): test-data coverage, A only, minimum (text 4.4)")
    check("N-383", max(e["tcr"] for e in g2), "6 (guard): test-data coverage, A only, maximum (text 4.4)")
    check("N-384", min(e["tcr"] for e in g4), "6 (guard): test-data coverage, neither, minimum (text 4.4)")
    check("N-385", max(e["tcr"] for e in g4), "6 (guard): test-data coverage, neither, maximum (text 4.4)")
    check("N-388", min(e["thr_r"] for e in g2), "6 (guard): recalculated threshold, A only, minimum (text 4.4)")
    check("N-389", max(e["thr_r"] for e in g2), "6 (guard): recalculated threshold, A only, maximum (text 4.4)")
    if any(e["thr_r"] != 0 for e in g4):
        stop("recalculated thresholds of the four-candidate cells are not all 0")
    check("N-390", sum(e["thr_r"] for e in g4), "6 (guard): recalculated threshold, neither (text 4.4)")
    check("N-386", min(e["thr"] for e in cells), "6 (guard): fixed thresholds, minimum (text 4.4)")
    check("N-387", max(e["thr"] for e in cells), "6 (guard): fixed thresholds, maximum (text 4.4)")
    check("N-022", sum(e["tcr"] < e["thr"] for e in cells),
          "6 lower tag line: cells whose test-data coverage under the real-data rule is below the fixed threshold "
          "(printed «22 of 24»)")
    check("N-392", sum(e["n_cand"] > 1 for e in cells), "6 (guard): cells with more than one candidate")
    if sorted(e["key"] for e in cells if e["tcr"] >= e["thr"]) != sorted(e["key"] for e in by_n[1]):
        stop("only the cells with both coils given should reach the fixed threshold in test data")
    odi_ref = by_n[1][1]
    same(PS.get("sym", "OXF1PRI~OXF2PRI|ODI", "threshold_true_combo_ours"), odi_ref["thr_r"], "symmetric threshold (d1_post_sym)")
    same(by_n[1][0]["thr_r"], by_n[1][0]["thr"], "NDI recalculated threshold equals the fixed one")
    check("N-257", None, "6 (guard): nine generator assumptions k (Section 3.4; no longer in the caption since fix 1)")
    coils = sorted({int(x) for e in cells for cmb in d1c[e["key"]]["candidates"] for x in cmb.split("|")})
    if coils != [32, 64]:
        stop(f"candidate coils {coils}")

    # ---------------- Fig. 6: Table S5 of the supplement (task: the tracks reproduce Table S5 for all 24 cells)
    s5 = parse_supp(SUPP_MD)
    seen, n_with = set(), 0
    for row in s5["rows"]:
        pair = tuple(row["Card pair"].split("~"))
        key = f"{row['Card pair']}|{row['Comp.']}"
        if key not in cj or len(pair) != 2:
            stop(f"Table S5 {row['Cell']}: {key} is not a calibration-test cell")
        seen.add(key)
        a, b = pair
        where = f"Table S5 {row['Cell']} ({key}, line {row['_line']})"
        if not printed_is(row["ΔTE (ms)"], acq[b]["te"] - acq[a]["te"], signed=True, where=f"{where} ΔTE"):
            stop(f"{where}: ΔTE {row['ΔTE (ms)']} is not TE({b}) − TE({a}) = {acq[b]['te'] - acq[a]['te']}")
        if not printed_is(row["ΔTR (s)"], acq[b]["tr"] - acq[a]["tr"], signed=True, where=f"{where} ΔTR"):
            stop(f"{where}: ΔTR {row['ΔTR (s)']} is not TR({b}) − TR({a}) = {acq[b]['tr'] - acq[a]['tr']}")
        if row["Same scheme"] != ("yes" if acq[a]["scheme"] == acq[b]["scheme"] else "no"):
            stop(f"{where}: «Same scheme» {row['Same scheme']} disagrees with the scheme track")
        both = acq[a]["given"] is not None and acq[b]["given"] is not None
        if row["Coils given"] != ("yes" if both else "no"):
            stop(f"{where}: «Coils given» {row['Coils given']} disagrees with the coil track")
        c = cj[key]
        if row["Comp."] == "NDI" and row["Role"] != "unblinded" or row["Comp."] == "ODI" and row["Role"] != "blinded":
            stop(f"{where}: role {row['Role']}")
        if c["verdict"] == "자료 없음":
            if (row["Coverage"], row["Threshold"], row["Result"], row["Subjects"]) != ("—", "—", "no data", "0"):
                stop(f"{where}: a cell without data is not printed as such")
            continue
        n_with += 1
        e = next(x for x in cells if x["key"] == key)
        if not printed_is(row["Coverage"], e["cov"], where=f"{where} Coverage") or \
                not printed_is(row["Threshold"], e["thr"], where=f"{where} Threshold"):
            stop(f"{where}: Coverage / Threshold {row['Coverage']} / {row['Threshold']} disagree with judgements.json")
        if int(row["Subjects"]) != e["n_subj"] or row["Result"] != ("pass" if e["passed"] else "fail"):
            stop(f"{where}: Subjects / Result disagree with judgements.json")
        br = s5["b"].get(row["Cell"])
        if br is None or br["pair_comp"] != f"{row['Card pair']} · {row['Comp.']}":
            stop(f"{where}: no matching row in Table S5 (b)")
        combos = sorted(f"{x}|{y}" for x, y in br["combos"])
        if combos != list(e["cands"]) or len(combos) != {"filled": 1}.get(acq[a]["chip"], 2) * \
                {"filled": 1}.get(acq[b]["chip"], 2):
            stop(f"{where}: Table S5 (b) combinations {combos} disagree with the coils given ({e['cands']})")
    if seen != set(cj) or n_with != len(cells) or len(s5["rows"]) != len(cj):
        stop(f"Table S5 has {len(s5['rows'])} rows for {len(cj)} cells ({n_with} with data)")
    s5_check = (f"Table S5 reproduced ({SUPP_MD.name}, sha256 {s5['sha'][:12]}): {len(s5['rows'])} rows "
                f"({n_with} with data, {len(s5['rows']) - n_with} without) — ΔTE, ΔTR, same scheme, coils given and the "
                "coil combinations of Table S5 (b) from the tracks; Coverage, Threshold, Subjects, Role, Result from "
                f"judgements.json; printed half to even at an exact tie: {len(S5_TIES)}"
                + (" (" + "; ".join(S5_TIES) + ")" if S5_TIES else ""))

    # ---------------- Fig. 7 (a) and (b): the Prisma pair
    meas = drawn("7a", "measured mean logit-NDI difference, 540 units (filled circle)", C1B, "sets", "all", "x", "mean")
    meas_ci = drawn("7a", "its subject-bootstrap 95 % interval (line)", C1B, "sets", "all", "x", "boot95")
    check("N-341", meas_ci[0], "7 (guard): lower end of the measured interval (drawn)")
    check("N-342", meas_ci[1], "7 (guard): upper end of the measured interval (drawn)")
    n_units = drawn("7a", "units of the measured difference (row label)", C1B, "sets", "all", "n_units")
    n_subj = drawn("7a", "participants of the measured difference (guard)", C1B, "sets", "all", "n_subjects")
    exp = {}
    for combo, lab in (("32|64", "recorded coils"), ("32|32", "32-channel setting at both")):
        exp[combo] = {m: drawn("7a", f"simulator expectation, {lab}, model {m} (gray tick)", C1, "expected", combo, "x", m)
                      for m in MODELS}
        same(C1.get("expected", combo, "x", "s3_clamped15"), exp[combo]["primary"], f"s3 = primary ({combo})")
        if max(exp[combo], key=exp[combo].get) != "primary":
            stop(f"the main model is not the largest expectation ({combo})")
    check("N-295", len(MODELS), "7 (guard): three regression models")
    cov_rec, cov_eq, thr_b = {}, {}, {}
    for comp in ("NDI", "ODI"):
        k = f"OXF1PRI~OXF2PRI|{comp}"
        cov_rec[comp] = drawn("7b", f"coverage with the recorded coils, {comp} (mark)", J, "V4", "cells", k, "coverage")
        same(SU.get("cells", f"32|64|{comp}", "coverage"), cov_rec[comp], f"suff_c3 32|64 {comp} = calibration-test coverage")
        thr_b[comp] = drawn("7b", f"fixed threshold, {comp} (line)", J, "V4", "cells", k, "threshold")
        for combo in ("32|32", "64|64"):
            cov_eq[(combo, comp)] = drawn("7b", f"coverage with equal settings {combo}, {comp} (mark)", SU, "cells",
                                          f"{combo}|{comp}", "coverage")
    same(thr_b["NDI"], thr_b["ODI"], "fixed thresholds of the two OXF1PRI~OXF2PRI cells")
    check("N-415", thr_b["ODI"], "7b: fixed threshold of both cells (printed)")
    check("N-423", min(cov_eq.values()), "7b: coverage with equal settings, minimum (printed)")
    check("N-424", max(cov_eq.values()), "7b: coverage with equal settings, maximum (printed)")
    if not all(v >= thr_b["NDI"] for v in cov_eq.values()) or any(v >= thr_b["NDI"] for v in cov_rec.values()):
        stop("fill of 7b: equal settings must reach the threshold and the recorded coils must not")
    thr_sym = drawn("7b", "ODI threshold of the equal-setting analysis (gray ticks)", PS, "sym", "OXF1PRI~OXF2PRI|ODI",
                    "threshold_true_combo_ours")
    check("N-416", thr_sym, "7b: ODI threshold of the equal-setting analysis (printed)")
    same(PS.get("sym", "OXF1PRI~OXF2PRI|NDI", "threshold_true_combo_ours"), thr_b["NDI"], "NDI symmetric = fixed threshold")
    for combo in ("32|32", "64|64"):
        same(SU.get("cells", f"{combo}|ODI", "threshold"), thr_sym, f"suff_c3 {combo} ODI threshold")
        same(SU.get("cells", f"{combo}|NDI", "threshold"), thr_b["NDI"], f"suff_c3 {combo} NDI threshold")
    m24 = re.search(r"대칭 문턱\(NDI ([\d.]+) · ODI ([\d.]+)\)", PL.text.splitlines()[1211])
    if not m24 or float(m24.group(1)) != thr_b["NDI"] or float(m24.group(2)) != thr_sym:
        stop("plan line 1212 does not fix NDI 0.80 · ODI 0.75 for the equal-setting analysis")
    for (combo, comp), cval in cov_eq.items():
        t_eq = thr_sym if comp == "ODI" else thr_b["NDI"]
        if (cval >= thr_b[comp]) != (cval >= t_eq):
            stop(f"fill of {combo} {comp} would change with the equal-setting threshold")
    if cov_rec["ODI"] >= thr_sym:
        stop("recorded ODI coverage reaches 0.75")

    # ---------------- Fig. 7 (c): the noise floor
    dev = {}
    for d in SIEMENS:
        e = dict(R=drawn("7c", f"{d}: session-median R (filled circle)", BS, "devices", d, "R_med"),
                 ci=drawn("7c", f"{d}: session-bootstrap 95 % interval of R (line)", BS, "devices", d, "R_ci"),
                 L=drawn("7c", f"{d}: effective L (printed beside the mark)", BS, "devices", d, "L_eff", "없음", "device", 0),
                 rice=KS.get("calibration", d, "없음", "1.0", "med"))        # Supplementary Table S10 (f), not drawn
        if not BS.get("devices", d, "separable") or not BA.get("devices", d, "separable"):
            stop(f"{d} should be determinable")
        same(BA.get("devices", d, "R_med"), e["R"], f"{d} R (all-scanner file)")
        for j in (0, 1):
            same(BA.get("devices", d, "R_ci")[j], e["ci"][j], f"{d} R interval (all-scanner file)")
        if not (e["ci"][0] <= e["R"] <= e["ci"][1]):
            stop(f"{d}: R outside its interval")
        dev[d] = e
    for i_, d in (("N-459", "OXF1PRI"), ("N-460", "OXF2PRI"), ("N-461", "OXF3TRI")):
        check(i_, dev[d]["rice"], f"7 (guard): {d} known-answer Rician R, 50 volumes (Supplementary Table S10 (f))")
    for d in OTHERS:
        e = dict(R=drawn("7c", f"{d}: session-median R (open circle)", BA, "devices", d, "R_med"),
                 ci=drawn("7c", f"{d}: session-bootstrap 95 % interval of R (line)", BA, "devices", d, "R_ci"))
        if BA.get("devices", d, "separable") or BA.get("devices", d, "vendor_measurable"):
            stop(f"{d} should be not determinable")
        dev[d] = e
    line = PL.text.splitlines()[1637]
    m1 = re.search(r"L ([\d.]+(?: · [\d.]+)+) 의 평균/SD", line)
    m2 = re.search(r"이론 ([\d.]+(?: · [\d.]+)+)", line)
    if not (m1 and m2):
        stop("plan line 1638 does not hold the chi reference values")
    l_list = [float(x) for x in m1.group(1).split(" · ")]
    theory = dict(zip(l_list, [float(x) for x in m2.group(1).split(" · ")]))
    if not {1.0, 4.0, 8.0} <= set(theory):
        stop(f"plan line 1638 lists L {l_list}, not 1, 4 and 8")
    ref = {L: theory[L] for L in (1.0, 4.0, 8.0)}
    for L, ids in ((1.0, ("N-462", "N-458")), (4.0, ("N-463",)), (8.0, ("N-464",))):
        for i in ids:
            check(i, ref[L], f"7c: large-sample mean/SD of a central chi floor, L {L:g} (dashed line)")
        drawn_md("7c", f"large-sample (population) mean/SD of a central chi floor, L {L:g} (dashed line)", PL, 1638, ref[L],
                 ids)
    mc = re.search(r"COIL_L = \{32: (\d+), 64: (\d+)\}", CF.text.splitlines()[50])
    if not mc:
        stop("s4cfg.py line 51 does not hold COIL_L")
    coil_l = {32: int(mc.group(1)), 64: int(mc.group(2))}
    check("N-429", coil_l[32], "7c: L of the 32-channel coil in the simulator (printed)")
    check("N-430", coil_l[64], "7c: L of the 64-channel coil in the simulator (printed)")
    check("N-276", [32, coil_l[32], 64, coil_l[64]], "7 (guard): simulator coil L")
    check("N-507", None, "7c: b = 2000 s/mm² volumes (printed as the shell name)")
    check("N-508", None, "7c: 50 volumes per session (printed)")
    for i in ("N-431", "N-432", "N-433", "N-434"):
        check(i, None, "7 (guard): effective-L ratios (text 4.4)")
    dec = BS.get("M_prime", "decision")
    if (dec["gate"], dec["side"], dec["size"]) != (True, "(가) 쪽", "크기 맞음"):
        stop(f"M′ decision {dec} is not «equal floors, expected size»")
    same(2 ** BS.get("M_prime", "contrasts", "c1", "c"), float(LEDGER["N-431"]["원 값"]), "N-431 = 2^c1")
    return dict(srcs_judge=srcs_judge, srcs7=srcs7, cells=cells, by_n=by_n, nodata=nodata, ans_ii=ans_ii, acq=acq,
                order=order, s5=s5, s5_check=s5_check, D1=D1, meas=meas, meas_ci=meas_ci, n_units=n_units, n_subj=n_subj,
                exp=exp, cov_rec=cov_rec, cov_eq=cov_eq, thr_b=thr_b["ODI"], thr_sym=thr_sym, dev=dev, ref=ref,
                coil_l=coil_l, PL=PL)


# ----------------------------------------------------------------------------- drawing helpers (shared with make_F7.py)
def mark(ax, x, y, comp=None, filled=True, color=INK, marker=None, ms=None, z=5, mew=0.8):
    ax.plot([x], [y], ls="none", marker=marker or MRK[comp], ms=ms or MS[comp], mfc=color if filled else SURF,
            mec=color, mew=mew, zorder=z, clip_on=False)


def vtick(ax, x, y, half, color, z=4, lw=1.0):
    ax.plot([x, x], [y - half, y + half], color=color, lw=lw, solid_capstyle="butt", zorder=z)


def width_mm(fig, t):
    fig.canvas.draw()
    bb = t.get_window_extent(fig.canvas.get_renderer())
    return bb.width / fig.dpi * 25.4


def tag(fig, cv, x, y, s, dashed):
    """Role tag: upper-case words in a box; solid border = test fixed before computing, dashed = after the verdict."""
    t = cv.text(x, y, s, fontsize=PT, color=INK, ha="left", va="center",
                bbox=dict(boxstyle="square,pad=0.28", fc=SURF, ec=INK2, lw=0.5, ls=TAG_DASHED if dashed else TAG_SOLID))
    return t


def header(fig, cv, x, y, letter, title, tag_text, dashed):
    """Panel letter, short title and role tag on one line; returns nothing."""
    cv.text(x, y, f"({letter})", fontsize=PTL, fontweight="bold", color=INK, ha="left", va="center")
    t = cv.text(x + 6.4, y, title, fontsize=PT, color=INK, ha="left", va="center")
    tag(fig, cv, x + 6.4 + width_mm(fig, t) + 2.6, y, tag_text, dashed)


def num(x):
    return f"{x:g}".replace("-", "−")


def signed(t):
    return "0" if abs(t) < 1e-12 else (f"{t:+.1f}".replace("-", "−"))


# ----------------------------------------------------------------------------- layout of Fig. 6 (mm; depth from the top)
CS = 6.5                                  # matrix cell = row pitch of the tracks
ROW0 = 19.0                               # top edge of the first row
MX = 104.0                                # left edge of the matrix
NAME_R = 19.0                             # acquisition names, right-aligned
VEND_C = 25.6
TE_X, TE_LIM, TE_TICKS = (32.5, 55.0), (62.0, 102.0), (70, 80, 90, 100)
TR_X, TR_LIM, TR_TICKS = (58.5, 78.5), (2.5, 11.0), (4, 6, 8, 10)
COIL_C, SCH_C = 85.2, 96.2
DNORM = TwoSlopeNorm(vcenter=0.0, vmin=-0.8, vmax=0.8)       # symmetric diverging scale (asserted to hold every value)
if DNORM.vmax != -DNORM.vmin:
    sys.exit("STOP F6/F7: the diverging scale must be symmetric about 0 (equal colour strength for equal distances)")
CB_TICKS = (-0.8, -0.4, 0.0, 0.4, 0.8)
CB = (161.5, 25.0, 2.4, 31.0)             # colour bar: x, top, width, height
CHIP_W, CHIP_H = 5.6, 3.5
CHECK = "✓"                          # ✓, set in DejaVu Sans (Liberation Sans has no check mark)


def oklab_L(rgba):
    return FC.oklab(matplotlib.colors.to_hex(rgba))[0]


def tick_color(fc):
    return "white" if oklab_L(fc) < 0.62 else INK


def draw6(v):
    order, acq = v["order"], v["acq"]
    n = len(order)
    bottom = ROW0 + n * CS
    y_coil, y_ndi = bottom + 16.0, bottom + 21.4          # fix 1: 0.8 mm lower, clear of the rotated column labels
    H = math.ceil(y_ndi + 3.2)
    fig = fs.new_fig(190, H)
    cv = fs.canvas(fig)
    Y = lambda d: H - d  # noqa: E731
    ymid = lambda i: Y(ROW0 + (i + 0.5) * CS)  # noqa: E731
    tx = lambda t: TE_X[0] + (t - TE_LIM[0]) / (TE_LIM[1] - TE_LIM[0]) * (TE_X[1] - TE_X[0])  # noqa: E731
    rx = lambda t: TR_X[0] + (t - TR_LIM[0]) / (TR_LIM[1] - TR_LIM[0]) * (TR_X[1] - TR_X[0])  # noqa: E731
    vals = [e["cov"] - e["thr"] for e in v["cells"]] + [e["tcr"] - e["thr"] for e in v["cells"]]
    if not all(DNORM.vmin <= x <= DNORM.vmax for x in vals):
        stop(f"a value outside the colour scale: {min(vals)} … {max(vals)}")
    if not all(TE_LIM[0] < acq[k]["te"] < TE_LIM[1] and TR_LIM[0] < acq[k]["tr"] < TR_LIM[1] for k in acq):
        stop("a TE or TR outside its track")

    # ---- role tags with the unfavourable results kept visible (S39b §2)
    t1 = tag(fig, cv, 2.0, Y(3.6), "BLINDED CONFIRMATION", dashed=False)
    t2 = tag(fig, cv, 2.0, Y(9.2), "AFTER THE VERDICT", dashed=True)
    xd = 2.0 + max(width_mm(fig, t1), width_mm(fig, t2)) + 3.2
    passed = show("6 upper tag line", ["N-266"], check("N-266"))
    cv.text(xd, Y(3.6), f"Upper triangle: real data. ODI {passed} passed; the NDI cell (unblinded check, below) failed.",
            fontsize=PT, color=INK, ha="left", va="center")
    below = show("6 lower tag line", ["N-022", "N-021"], f"{check('N-022')} of {check('N-021')}")
    cv.text(xd, Y(9.2), f"Lower triangle: the judge's own simulated test data, scored with the real-data rule. {below} cells "
            "below the fixed threshold.", fontsize=PT, color=INK, ha="left", va="center")

    # ---- track headers
    yh = Y(ROW0 - 1.0)
    for x, s in ((VEND_C, "vendor"), ((TE_X[0] + TE_X[1]) / 2, "TE (ms)"), ((TR_X[0] + TR_X[1]) / 2, "TR (s)"),
                 (COIL_C, "receive\ncoil"), (SCH_C, "scheme")):
        cv.text(x, yh, s, fontsize=PT, color=INK2, ha="center", va="bottom", linespacing=1.0)

    # ---- tracks
    yb = Y(bottom + 0.7)
    for (x0, x1), f, ticks in ((TE_X, tx, TE_TICKS), (TR_X, rx, TR_TICKS)):
        for t in ticks:
            cv.plot([f(t), f(t)], [yb, Y(ROW0)], color=GRID, lw=0.5, zorder=1)
            cv.plot([f(t), f(t)], [yb, yb - 0.8], color=AXIS, lw=0.5, zorder=1)
            cv.text(f(t), yb - 1.2, f"{t:g}", fontsize=PT, color=INK2, ha="center", va="top")
        cv.plot([x0, x1], [yb, yb], color=AXIS, lw=0.5, solid_capstyle="butt", zorder=1)
    for i, nm in enumerate(order):
        a, y = acq[nm], ymid(i)
        cv.text(NAME_R, y, nm, fontsize=PT, color=INK, ha="right", va="center")
        cv.text(VEND_C, y, a["vendor"], fontsize=PT, color=INK2, ha="center", va="center")
        for (x0, x1), xx in ((TE_X, tx(a["te"])), (TR_X, rx(a["tr"]))):
            cv.plot([x0, x1], [y, y], color=FC.GRID, lw=0.5, solid_capstyle="butt", zorder=1.5)
            cv.plot([xx], [y], ls="none", marker="o", ms=3.6, mfc=INK2, mec="white", mew=0.5, zorder=3)
        if a["chip"] == "none":
            cv.text(COIL_C, y, "–", fontsize=PT, color=MUTED, ha="center", va="center")
        else:
            filled = a["chip"] == "filled"
            cv.add_patch(FancyBboxPatch((COIL_C - CHIP_W / 2, y - CHIP_H / 2), CHIP_W, CHIP_H,
                                        boxstyle="round,pad=0,rounding_size=0.9", fc=INK2 if filled else SURF, ec=INK2,
                                        lw=0.6, zorder=2))
            cv.text(COIL_C, y, str(a["given"] if filled else a["rec"]), fontsize=PT, color=SURF if filled else INK2,
                    ha="center", va="center", zorder=3)
        cv.text(SCH_C, y, a["scheme"], fontsize=PT, color=INK2, ha="center", va="center")

    # ---- matrix
    drawn_at = {}

    def cell(r, c, value=None, thr=None, hatch=False, outline=False, x0=None, ytop=None, size=CS):
        x0 = MX + c * CS if x0 is None else x0
        ytop = ROW0 + r * CS if ytop is None else ytop
        if hatch:
            cv.add_patch(Rectangle((x0 + 0.35, Y(ytop + size) + 0.35), size - 0.7, size - 0.7, fc=SURF, ec=FC.MUTED,
                                   lw=0.5, hatch="////", zorder=2))
            fc = None
        else:
            fc = FC.DIV(DNORM(value - thr))
            cv.add_patch(Rectangle((x0, Y(ytop + size)), size, size, fc=fc, ec=SURF, lw=1.2, zorder=2))
        if outline:
            cv.add_patch(Rectangle((x0 + 0.55, Y(ytop + size) + 0.55), size - 1.1, size - 1.1, fc="none", ec=INK, lw=1.2,
                                   zorder=4))
        if not hatch and value >= thr:
            cv.text(x0 + size / 2, Y(ytop + size / 2), CHECK, fontfamily="DejaVu Sans", fontsize=8.5,
                    color=tick_color(fc), ha="center", va="center", zorder=5)
        return fc

    for i in range(n):
        cv.add_patch(Rectangle((MX + i * CS + 0.35, Y(ROW0 + (i + 1) * CS) + 0.35), CS - 0.7, CS - 0.7, fc=SURF,
                               ec=GRID, lw=0.5, zorder=2))
    n_tick = 0
    for e in v["cells"]:
        if e["comp"] != "ODI":
            continue
        i, j = sorted((order.index(e["pair"][0]), order.index(e["pair"][1])))
        pr = e["n_cand"] == 1
        cell(i, j, e["cov"], e["thr"], outline=pr)
        cell(j, i, e["tcr"], e["thr"], outline=pr)
        n_tick += (e["cov"] >= e["thr"]) + (e["tcr"] >= e["thr"])
        for rc in ((i, j), (j, i)):
            drawn_at[rc] = drawn_at.get(rc, 0) + 1
    for p in v["nodata"]:
        i, j = sorted((order.index(p[0]), order.index(p[1])))
        cell(i, j, hatch=True)
        cell(j, i, hatch=True)
        for rc in ((i, j), (j, i)):
            drawn_at[rc] = drawn_at.get(rc, 0) + 1
    if sorted(drawn_at) != sorted((r, c) for r in range(n) for c in range(n) if r != c) or set(drawn_at.values()) != {1}:
        stop("the matrix cells are not drawn exactly once each")
    if n_tick != sum(e["tcr"] >= e["thr"] for e in v["cells"] if e["comp"] == "ODI") or \
            any(e["cov"] >= e["thr"] for e in v["cells"]):
        stop("check marks do not follow the values")
    labels = []
    for j, nm in enumerate(order):
        labels.append(cv.text(MX + (j + 0.5) * CS + 1.1, Y(bottom + 0.9), nm, fontsize=PT, color=INK, rotation=50,
                              rotation_mode="anchor", ha="right", va="top"))

    # ---- colour bar and keys (right)
    cx, ctop, cw, ch = CB
    cv.text(cx - 2.6, Y(ctop - 1.6), "Coverage of zero\nminus the cell's\nfixed threshold", fontsize=PT, color=INK,
            ha="left", va="bottom", linespacing=1.05)
    cax = fs.ax_mm(fig, cx, H - ctop - ch, cw, ch)
    cb = fig.colorbar(ScalarMappable(norm=DNORM, cmap=FC.DIV), cax=cax)
    cb.outline.set_linewidth(0.5)
    cb.outline.set_edgecolor(AXIS)
    cax.yaxis.set_major_locator(FixedLocator(CB_TICKS))
    cax.yaxis.set_minor_locator(NullLocator())
    cax.yaxis.set_major_formatter(FuncFormatter(lambda t, _: signed(t)))
    cax.tick_params(axis="y", labelsize=PT, length=2.0, width=0.5, colors=INK2, pad=1.5)
    ky = ctop + ch + 5.0
    kx = cx - 2.6
    ks = 4.2
    pr_odi = v["by_n"][1][1]                         # key swatch: the colour of the one ✓ cell of the matrix (stored values)
    if pr_odi["comp"] != "ODI" or pr_odi["tcr"] < pr_odi["thr"]:
        stop("the key swatch must be the Prisma ODI cell at or above its threshold")
    cell(0, 0, pr_odi["tcr"], pr_odi["thr"], x0=kx, ytop=ky - ks / 2, size=ks)
    cv.text(kx + ks + 1.6, Y(ky), "at or above the\nfixed threshold", fontsize=PT, color=INK, ha="left", va="center",
            linespacing=1.05)
    cv.add_patch(Rectangle((kx, Y(ky + 7.0 + ks / 2)), ks, ks, fc=FC.NEUTRAL, ec=SURF, lw=1.0, zorder=2))
    cv.add_patch(Rectangle((kx + 0.45, Y(ky + 7.0 + ks / 2) + 0.45), ks - 0.9, ks - 0.9, fc="none", ec=INK, lw=1.2,
                           zorder=4))
    cv.text(kx + ks + 1.6, Y(ky + 7.0), "the Prisma pair", fontsize=PT, color=INK, ha="left", va="center")
    cell(0, 0, hatch=True, x0=kx, ytop=ky + 12.4 - ks / 2, size=ks)
    cv.text(kx + ks + 1.6, Y(ky + 12.4), "no data", fontsize=PT, color=INK, ha="left", va="center")

    # ---- receive-coil key and the NDI cell (footnote lines)
    y = Y(y_coil)
    t = cv.text(2.0, y, "Receive coil:", fontsize=PT, color=INK, ha="left", va="center")
    x = 2.0 + width_mm(fig, t) + 2.0
    for kind, words in (("filled", "given to the judge"), ("open", "recorded but not given"), ("none", "not recorded")):
        if kind == "none":
            t = cv.text(x + 1.0, y, "–", fontsize=PT, color=MUTED, ha="center", va="center")
            x += 2.6
        else:
            cv.add_patch(FancyBboxPatch((x, y - CHIP_H / 2), CHIP_W, CHIP_H, boxstyle="round,pad=0,rounding_size=0.9",
                                        fc=INK2 if kind == "filled" else SURF, ec=INK2, lw=0.6, zorder=2))
            x += CHIP_W + 1.4
        t = cv.text(x, y, words, fontsize=PT, color=INK, ha="left", va="center")
        x += width_mm(fig, t) + 4.5
    # the open chip and the dash come from the sidecar records (rv_4a_coil.json), read after the verdict: graded here
    cv.text(x - 1.5, y, "(sidecar records, read after the verdict)", fontsize=PT, color=INK, ha="left", va="center")
    ndi = v["by_n"][1][0]
    if ndi["comp"] != "NDI":
        stop("the NDI cell is not first in the group with both coils given")
    y = Y(y_ndi)
    t = cv.text(2.0, y, "NDI cell, OXF1PRI~OXF2PRI (unblinded check; not in the matrix):", fontsize=PT, color=INK,
                ha="left", va="center")
    x = 2.0 + width_mm(fig, t) + 2.4
    ms_ = 4.2
    cell(0, 0, ndi["cov"], ndi["thr"], outline=True, x0=x, ytop=y_ndi - ms_ / 2, size=ms_)
    t = cv.text(x + ms_ + 1.4, y, "real data", fontsize=PT, color=INK, ha="left", va="center")
    x += ms_ + 1.4 + width_mm(fig, t) + 4.0
    cell(0, 0, ndi["tcr"], ndi["thr"], outline=True, x0=x, ytop=y_ndi - ms_ / 2, size=ms_)
    cv.text(x + ms_ + 1.4, y, "test data, after the verdict (Section 4.4; Supplementary Table S5)", fontsize=PT,
            color=INK, ha="left", va="center")
    v["_labels"] = labels
    return fig


def text_overlaps(fig):
    """Pairs of printed texts whose boxes overlap by more than 0.25 mm both ways (none allowed). Texts rotated by the same
    non-zero angle (the column labels) are compared by their perpendicular spacing instead (checked in parallel_ok)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    px_mm = fig.dpi / 25.4
    boxes = []
    for t in fig.findobj(mtext.Text):
        if not t.get_visible() or not t.get_text().strip():
            continue
        try:
            bb = t.get_window_extent(r)
        except RuntimeError:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        boxes.append((t.get_text(), bb, round(t.get_rotation(), 6)))
    tol = 0.25 * px_mm
    bad = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            if boxes[i][2] == boxes[j][2] != 0:
                continue
            a, b = boxes[i][1], boxes[j][1]
            ox = min(a.x1, b.x1) - max(a.x0, b.x0)
            oy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ox > tol and oy > tol:
                bad.append((boxes[i][0], boxes[j][0], round(ox / px_mm, 2), round(oy / px_mm, 2)))
    return bad


def parallel_ok(fig, labels, pitch_mm):
    """Parallel rotated labels at a fixed pitch do not touch when the pitch across their direction exceeds their height."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    for t in labels:
        rot = math.radians(t.get_rotation())
        w, h, d = r.get_text_width_height_descent(t.get_text(), t.get_fontproperties(), ismath=False)
        if pitch_mm * math.sin(rot) <= (h + d) / fig.dpi * 25.4 + 0.25:
            stop(f"rotated labels at pitch {pitch_mm} mm would touch («{t.get_text()}»)")


# ----------------------------------------------------------------------------- caption
SRC6 = ("[src: `decomp/results/remeasure/judgements.json` `V4.cells` · `V4.summary` · ledger copies of "
        "`remeasure/explore/d1/d1_rulethr.json` `cells` (`candidates` · `threshold_stored` · `test_coverage_true_combo_kmed` · "
        "`test_coverage_rule_kmed` · `n_clusters`), `remeasure/explore/d1/d1_post_sym.json` `sym`, "
        "`remeasure/explore/rv/rv_4a_coil.json` `by_device` · `sessions`, `decomp/results/stage3/cards/cards.json` (`dTE` · "
        "`dTR_ms` · `design_differs`), `decomp/s4/plan.py:17–22` `ONH` and `decomp/stage1/analyze3.py` `NOT2ING_FIXED` · "
        "drawn by `decomp/manuscript/figures/src/make_F6.py`, which stops unless every printed number equals its ledger row, "
        "every drawn value equals its field in the ledger copies and, where one exists, its ledger row, and the tracks "
        "reproduce ΔTE, ΔTR, scheme and coil combinations of Supplementary Table S5]{.src}")
# Fix 1 (2026-11-06): one table reference for both S5 variants (direction S42 §1.3 asks for columns; the default writes
# part (e)), so the caption stays right whichever variant is adopted; «Supplementary Table S5» holds part (e) too.
S5_REF = "Supplementary Table S5"
# Fix 1: the outline names the pair (two cells: ODI in the matrix, NDI below it) in the body's own words (Section 4.6:
# «Only in OXF1PRI~OXF2PRI, with both coils given to the judge, were the thresholds set under the rule applied to the
# data.»); the NDI line is graded (unblinded check); the scheme letters and the after-the-verdict coil records are
# defined. «median over the nine generator assumptions» is left to Section 3.4, which defines cell coverage so.
CAPTION = ("![**Fig. 6. The frozen judge on traveling heads.** Acquisitions ordered by TE, with vendor, TE, TR, receive coil "
           "(open chip and dash, sidecar records read after the verdict) and direction scheme (same letter, same scheme). "
           "Cells: coverage of zero minus the cell's fixed threshold, {N-380} ODI cells; upper triangle, real data "
           "(blinded); lower, after the verdict, the judge's own simulated test data scored with the real-data rule ("
           + S5_REF + "). Check mark, threshold met. Outlined, the Prisma pair, the only pair whose thresholds were set "
           "under the rule applied to the data; below the matrix, its NDI cell (unblinded check). Hatched, no data. "
           + SRC6 + "](figures/F6_judge_matrix.pdf){#fig:fig6 .wide width=100%}")


def caption(template):
    def sub(m):
        i = m.group(1)
        check(i, None, "caption")
        return f"{LEDGER[i]['원고 표기']} [{i}]"
    return re.sub(r"\{(N-\d{3})\}", sub, template)


def caption_words(cap):
    """Words of the caption without the src block and the ledger tags (the reviewer's count)."""
    t = re.sub(r"\[src:[^\]]*\]\{\.src\}", "", cap)
    t = re.sub(r"\[N-\d{3}\]", "", t)
    t = re.sub(r"\]\(figures/[^)]*\)\{[^}]*\}$", "", t.strip())
    return len(t.replace("![", "").split())


def ledger_check(path):
    """Run the project's ledger check on the caption with every region strict; stop unless it passes."""
    r = subprocess.run([sys.executable, "-I", "-B", str(LEDGER_DIR / "ledger_check.py"), str(path), "--strict-all", "1"],
                       capture_output=True, text=True)
    out = r.stdout.strip()
    if r.returncode != 0:
        stop(f"ledger_check on {path.name}: {out} {r.stderr.strip()}")
    return out.splitlines()[0] if out else r.stderr.strip()


# ----------------------------------------------------------------------------- Supplementary Table S5: after the verdict
TS5_SRC = ("[src: ledger copy of `remeasure/explore/d1/d1_rulethr.json` (sha256 {sha}…) `cells.<cell>` `candidates` · "
           "`test_coverage_rule_kmed` · `threshold_rule_ours` · `threshold_stored` (= `judgements.json` `V4.cells.<cell>.threshold`), "
           "rounded half up to two decimals and written by `decomp/manuscript/figures/src/make_F6.py`, which stops unless "
           "the minima and maxima per coil group equal the ledger rows of Section 4.4]{{.src}}")
TS5_MARK = ("[A 고침 · 2026-11-06 · v12 — 방향 S42 §1.3: 판정 뒤 칸별 값(시험 자료 덮음 · 다시 낸 문턱)을 표 S5 에 · "
            "값은 make_F6.py 가 장부 사본에서 씀(손으로 옮기지 않음)]")
# Fix 1 (2026-11-06, verifier finding 5): why the Prisma ODI row shows a recalculated threshold that differs from the fixed
# one although the rule is the same. In d1_c3.py mode_rulethr, a cell with one candidate combination gets the same
# coverage, «ok» and flag arrays under the real-data rule as under the true-combination rule, and the same random key,
# so threshold_rule_ours equals threshold_true_combo_ours (ledger N-416 holds the latter for this cell); the difference
# from threshold_stored therefore lies in the test cases of the recalculation. No number is printed in the sentence.
TS5_PRISMA = ("For OXF1PRI~OXF2PRI (one candidate combination), the rule applied to real data is the one its fixed "
              "thresholds were set under; a difference between its recalculated and fixed thresholds therefore comes from "
              "the test cases used for the recalculation, not from the rule.")
TS5_PRISMA_MARK = ("[A 제안 · 2026-11-06 · v12 그림 6 고침 1 — Prisma 쌍의 다시 낸 문턱이 고정 문턱과 다른 까닭: 후보 하나인 칸은 "
                   "실측 규칙과 참 조합 규칙이 같은 배열 · 같은 난수 열쇠(d1_c3.py mode_rulethr) · 같은 사례의 참 조합 문턱 "
                   "= 다시 낸 문턱(장부 N-416) · 저장 문턱과의 차이는 사례에서 · A 확인 뒤 이 표지 지움]")


def s5_values(v):
    """Per cell (Table S5 order): the two after-the-verdict values as printed, with their guards."""
    D1 = v["D1"]
    out = {}
    for row in v["s5"]["rows"]:
        key = f"{row['Card pair']}|{row['Comp.']}"
        e = next((x for x in v["cells"] if x["key"] == key), None)
        if e is None:
            out[row["Cell"]] = None
            continue
        t_c, t_r, t_f = (str(half_up(e["tcr"], 2)), str(half_up(e["thr_r"], 2)), str(half_up(e["thr"], 2)))
        el = f"Table S5 {row['Cell']} {key}"
        for what, x, t in (("test-data coverage", e["tcr"], t_c), ("recalculated threshold", e["thr_r"], t_r)):
            if t != f"{x:.2f}":
                S5_NEW_TIES.append(f"{row['Cell']} {what} {x!r} -> {t} (half up, the ledger's rule; Python's format "
                                   f"would print {x:.2f})")
        tabled(f"{el}: test-data coverage under the real-data rule, median over k", D1,
               ("cells", key, "test_coverage_rule_kmed"), e["tcr"], t_c)
        tabled(f"{el}: threshold recalculated with the real-data rule", D1, ("cells", key, "threshold_rule_ours"),
               e["thr_r"], t_r)
        if t_f != row["Threshold"]:
            stop(f"{el}: fixed threshold {t_f} is not Table S5's {row['Threshold']}")
        out[row["Cell"]] = dict(e=e, tcr=t_c, thr_r=t_r, thr=t_f)
    # the printed extremes per coil group are the ledger rows of Section 4.4 (two decimals)
    for n, (ilo, ihi, jlo, jhi) in ((2, ("N-382", "N-383", "N-388", "N-389")), (4, ("N-384", "N-385", "N-390", "N-390"))):
        g = [x for x in out.values() if x and x["e"]["n_cand"] == n]
        for ids, fld in (((ilo, ihi), "tcr"), ((jlo, jhi), "thr_r")):
            lo = min(g, key=lambda x: Decimal(x[fld]))[fld]
            hi = max(g, key=lambda x: Decimal(x[fld]))[fld]
            for want, got in ((LEDGER[ids[0]]["원고 표기"], lo), (LEDGER[ids[1]]["원고 표기"], hi)):
                if Decimal(want) != Decimal(got):
                    stop(f"Table S5 after the verdict, {n} candidates, {fld}: printed {got} but Section 4.6 has {want}")
    return out


def s5_edits(v, mode, vals):
    s5 = v["s5"]
    L = s5["lines"]
    src = TS5_SRC.format(sha=v["D1"].sha[:12])
    edits = []
    if mode == "columns":
        old = "\n".join(L[s5["head"]:s5["last"] + 1])
        new = [L[s5["head"]] + " Test data, real-data rule | Threshold, real-data rule |", L[s5["head"] + 1] + "------|------|"]
        for row in s5["rows"]:
            x = vals[row["Cell"]]
            new.append(row["_text"] + (f" {x['tcr']} | {x['thr_r']} |" if x else " — | — |"))
        edits.append(dict(id="F6-S5-1", target="supplement", old=old, new="\n".join(new),
                          note="Table S5: two after-the-verdict columns (test-data coverage under the real-data rule, median "
                               "over k; threshold recalculated with that rule), values written by make_F6.py"))
        anchor = "· coil `decomp/s4/plan.py:17–22`]{.src}"
        edits.append(dict(id="F6-S5-2", target="supplement", old=anchor, new=anchor + (
            " [After the verdict (post hoc), the last two columns give the coverage of zero of the judge's own simulated "
            "test data of the cell's acquisition type, restricted to cases whose true receive-coil combination is among the "
            "cell's candidates (Table S5 (b)) and scored with the rule applied to real data, as the median over the nine "
            "generator assumptions k, and the threshold recalculated from the same cases with that rule; both to two "
            "decimals, as in Section 4.6. They change no verdict (Section 4.6; Fig. 6, lower triangle). " + TS5_PRISMA +
            "]{.aedit} " + TS5_MARK + " " + TS5_PRISMA_MARK + " " + src),
            note="Table S5 notes: what the two new columns are"))
    elif mode == "part_e":
        rows = ["| Cell | Card pair · component | Candidate coil combinations | Test-data coverage, real-data rule | "
                "Threshold recalculated with that rule | Fixed threshold |",
                "|------|------------------|------|------|------|------|"]
        for row in s5["rows"]:
            x = vals[row["Cell"]]
            if x:
                rows.append(f"| {row['Cell']} | {row['Card pair']} · {row['Comp.']} | {x['e']['n_cand']} | {x['tcr']} | "
                            f"{x['thr_r']} | {x['thr']} |")
        block = ("Table: **Table S5 (e).** [After the verdict (post hoc): the judge's own simulated test data scored with "
                 "the rule applied to real data, for the 24 cells with data in the order of Table S5. Test-data coverage: "
                 "the coverage of zero in the test data of the cell's acquisition type, restricted to cases whose true "
                 "receive-coil combination is among the cell's candidate combinations (Table S5 (b)), a unit counting as "
                 "covered only if every candidate combination covered zero, as the median over the nine generator "
                 "assumptions k. Recalculated threshold: the 5th percentile of that coverage, from the same cases. Fixed "
                 "threshold: as in Table S5. Values to two decimals, as in Section 4.6; they change no verdict (Fig. 6, "
                 "lower triangle). " + TS5_PRISMA + "]{.aedit} " + TS5_MARK + " " + TS5_PRISMA_MARK + " " + src + "\n\n" +
                 "\n".join(rows) + "\n\n")
        anchor = "### Table S6 {.unnumbered}"
        edits.append(dict(id="F6-S5e-1", target="supplement", old=anchor, new=block + anchor,
                          note="new Table S5 (e) after Table S5 (d), same cell order as Table S5; values written by "
                               "make_F6.py"))
        anchor2 = "· coil `decomp/s4/plan.py:17–22`]{.src}"
        edits.append(dict(id="F6-S5e-2", target="supplement", old=anchor2, new=anchor2 + (
            " [After the verdict, Table S5 (e) gives for each cell the coverage of the judge's own simulated test data "
            "under the rule applied to real data and the threshold recalculated with that rule (Section 4.6; Fig. 6, lower "
            "triangle).]{.aedit} " + TS5_MARK),
            note="Table S5 notes: pointer to part (e)"))
    else:
        stop(f"S5_MODE {mode!r}")
    for e in edits:
        k = s5["text"].count(e["old"])
        if k != 1:
            stop(f"{e['id']}: «old» occurs {k} times in {SUPP_MD.name}")
    return edits


# ----------------------------------------------------------------------------- values file
def short(val):
    if isinstance(val, float):
        return repr(val)
    if isinstance(val, list):
        return "[" + ", ".join(short(x) for x in val) + "]"
    return repr(val)


def cell_md(s):
    """A Markdown table cell: escape the pipes that ledger keys carry (32|64, OXF1PRI~OXF2PRI|NDI)."""
    return str(s).replace("|", "\\|")


def values_md(fig_no, name, srcs, cap, check_line, audit_line, n_words, H, notes, extra=()):
    pre = str(fig_no)
    shown = [s for s in SHOWN if s["where"].startswith(pre)]
    drawn_ = [d for d in DRAWN if d["panel"].startswith(pre)]
    used = {i: [r for r in roles if r.startswith(pre) or r == "caption"] for i, roles in USED.items()}
    used = {i: r for i, r in used.items() if r}
    L = [f"# Figure {fig_no} — drawn values and visible numbers", "",
         f"Made by `make_F{fig_no}.py` (v12 revision of 2026-11-06; direction S42). Ledger: `drafts/x15/ledger/numbers.csv` "
         f"({len(LEDGER)} rows). The script stops unless every drawn value that the ledger also holds (same file, same "
         "sha256, same field) equals the row's raw value (원 값), every ledger row used re-derives from its source copy with "
         "`ledger_check.py`'s own functions and its 원고 표기 passes its rounding rule, and every printed number is the "
         "원고 표기 of a row whose raw value equals the drawn value (for a range: the drawn minimum and maximum; for a count: "
         "counted from the cells).",
         "", f"Audit: `{audit_line}` · figure 190 × {H} mm · caption {n_words} words without the src block and the tags.", "",
         "## Sources read", "", "| file (ledger name) | read from | sha256 |", "|---|---|---|"]
    for s in srcs.values():
        L.append(f"| `{s.rel}` | ledger copy `sources/{s.path.relative_to(SRC)}` | `{s.sha}` |")
    L += ["", "The remeasure files behind these ledger copies and the repository's `judgements.json` and `cards.json` were "
          "logged as re-read (drafts/x14/REREAD_LOG.jsonl) before the run.", "", "## Visible numbers in the figure", "",
          "Axis and colour-bar tick labels are exempt. Scanner and acquisition names (OXF1PRI, NOT2ING-89, …), vendor names, "
          "the scheme letters A and B, the coil names 32 and 64 in the receive-coil chips (32-channel and 64-channel "
          "coils) and section or table numbers are names, not values (exempt as in `ledger_check.py` EXCEPTIONS).", "",
          "| where | text | ledger ID(s) | 원고 표기 | raw value(s) |", "|---|---|---|---|---|"]
    for s in shown:
        L.append(f"| {cell_md(s['where'])} | {cell_md(s['text'])} | {' · '.join(s['ids'])} | "
                 f"{cell_md(' · '.join(LEDGER[i]['원고 표기'] for i in s['ids']))} | "
                 f"{cell_md(' · '.join(LEDGER[i]['원 값'] for i in s['ids']))} |")
    L += ["", "Notes:"] + [f"- {n}" for n in notes]
    L += ["", "## Drawn values", "",
          "Every value the script draws, with its source field. «no ledger row (drawn only)» marks a drawn value that the "
          "ledger does not hold under the same file and field; none of them is printed in the figure.", "",
          "| panel | element | file (sha256[:12]) | field | value | ledger ID |", "|---|---|---|---|---|---|"]
    n_new = 0
    for d in drawn_:
        lid = " · ".join(d["ids"]) if d["ids"] else "no ledger row (drawn only)"
        n_new += not d["ids"]
        L.append(f"| {d['panel']} | {cell_md(d['element'])} | `{d['file']}` ({d['sha']}) | `{cell_md(d['field'])}` | "
                 f"{cell_md(short(d['value']))} | {lid} |")
    L += ["", f"Drawn values: {len(drawn_)}; with a ledger row under the same field: {len(drawn_) - n_new}; "
          f"without a ledger row (drawn only, not printed): {n_new}.", ""]
    L += list(extra)
    L += ["## Ledger rows used (re-derived from their source copies)", "",
          "| ledger ID | 원고 표기 | raw value | field | used for |", "|---|---|---|---|---|"]
    for i in sorted(used):
        r = LEDGER[i]
        L.append(f"| {i} | {cell_md(r['원고 표기'])} | {cell_md(r['원 값'][:60])} | `{r['근거 파일']}` `{cell_md(r['근거 필드'][:90])}` | "
                 f"{cell_md('; '.join(sorted(set(used[i]))))} |")
    L += ["", f"## Numbers in the caption (F{fig_no}_caption.md)", "", "| text | ledger ID | 확인 상태 |", "|---|---|---|"]
    for m in re.finditer(r"\[(N-\d{3})\]", cap):
        L.append(f"| {cell_md(LEDGER[m.group(1)]['원고 표기'])} | {m.group(1)} | {LEDGER[m.group(1)]['확인 상태']} |")
    L += ["", "## Ledger check of the caption", "",
          f"`ledger_check.py F{fig_no}_caption.md --strict-all 1` (run by the script, which stops unless it exits 0):", "",
          f"    {check_line}", ""]
    return "\n".join(L) + "\n"


class Tee:
    def __init__(self):
        self.lines = []

    def write(self, s):
        sys.__stdout__.write(s)
        self.lines.append(s)

    def flush(self):
        sys.__stdout__.flush()


def finish(fig, fig_no, name, template, srcs, notes, extra=()):
    """Overlap check, caption with ledger check, save with the figstyle audit, values file; shared by make_F7.py.
    Returns the audit line."""
    cap = caption(template)
    n_words = caption_words(cap)
    bad = text_overlaps(fig)
    if bad:
        stop("overlapping texts: " + "; ".join(f"«{a[:30]}» × «{b[:30]}» ({ox} × {oy} mm)" for a, b, ox, oy in bad[:12]))
    if fs.OUTDIR != TARGET:
        stop(f"figstyle OUTDIR changed to {fs.OUTDIR}")
    tmp = OUT / f"F{fig_no}_caption.md.tmp"              # checked first, so a failing caption never replaces the last one
    tmp.write_text(cap + "\n", encoding="utf-8")
    try:
        check_line = ledger_check(tmp)
    except SystemExit:
        tmp.unlink(missing_ok=True)
        raise
    tmp.replace(OUT / f"F{fig_no}_caption.md")
    print("ledger check (caption):", check_line)
    print("caption words (no src, no tags):", n_words)
    H = round(fig.get_size_inches()[1] * 25.4)
    tee = Tee()
    sys.stdout = tee
    try:
        fs.save(fig, name)                               # refuses to write unless the audit passes
    finally:
        sys.stdout = sys.__stdout__
    audit_line = next((s.strip() for s in "".join(tee.lines).splitlines() if s.startswith("AUDIT")), "")
    (OUT / f"F{fig_no}_values.md").write_text(values_md(fig_no, name, srcs, cap, check_line, audit_line, n_words, H, notes,
                                                        extra), encoding="utf-8")
    print(f"printed numbers {len([s for s in SHOWN if s['where'].startswith(str(fig_no))])} · drawn values "
          f"{len([d for d in DRAWN if d['panel'].startswith(str(fig_no))])}")
    for ext in ("pdf", "png"):
        f = OUT / f"{name}.{ext}"
        print(f"sha256 {f.name} {hashlib.sha256(f.read_bytes()).hexdigest()}")
    return audit_line


def notes6(v):
    def rng(vals):
        return f"{min(vals):g} – {max(vals):g} ({len(vals)} distinct values)"
    a = v["ans_ii"]
    return [
        "Design: review direction S42 §1.3 (the prototype fig_v3_matrices.py is the layout reference; its values are not "
        "used). One matrix replaces v11 Fig. 5 (a) and (b): upper triangle = old (a) (real data), lower triangle = old (b) "
        "(test data under the real-data rule). The range over k (old (a) line) and the recalculated thresholds (old (b) "
        "gray ticks) are not drawn: the range stays guarded, the recalculated thresholds are printed per cell in "
        "Supplementary Table S5 (written by this script).",
        "Colour: figcolors.DIV on a symmetric scale, TwoSlopeNorm(vmin −0.8, vcenter 0, vmax +0.8) (symmetry asserted); "
        "every drawn value lies inside it (asserted). Fix 1 (2026-11-06, verifier finding 4): the prototype's scale "
        "(−0.8 … 0 … +0.2) gave each arm the full colour range, so the one favourable ODI cell (Prisma pair, lower "
        "triangle) was painted as strongly as a cell about four times as far below its threshold; on the symmetric scale "
        "equal distances from the threshold get equal colour strength, and the cells above the threshold are pale blue "
        "with ✓. ✓ marks exactly the values at or above the fixed threshold (asserted): none in the upper triangle, only "
        "the Prisma pair in the lower one (and the NDI cell in the footnote).",
        "Outline (fix 1, verifier finding 1): the Prisma pair has two cells, ODI (in the matrix, both triangles) and NDI "
        "(footnote line), both outlined; the caption now says «the only pair whose thresholds were set under the rule "
        "applied to the data», the body's own words (Section 4.6, [자료]), and grades the NDI line «(unblinded check)» "
        "(finding 2). Scheme letters are defined in the caption, «same letter, same scheme», as in Fig. 5 (finding 3).",
        "Supplementary Table S5 (fix 1, verifier finding 6): direction S42 §1.3 asks for the after-the-verdict values «as "
        "Table S5 columns»; the default (S5_MODE part_e) writes them as a new part (e) instead, and both variants are "
        "written (F6_supp_edits_columns.json, F6_supp_edits_part_e.json). Deviation from the direction, left to A. The "
        "caption cites «Supplementary Table S5» in both variants, so it stays right whichever is adopted. Both variants "
        "add one sentence on the Prisma pair (finding 5): with one candidate combination the real-data rule is the rule "
        "its fixed thresholds were set under, so the recalculated ODI threshold differs from the fixed one through the "
        "test cases, not the rule (premise asserted: the true-combination threshold of d1_post_sym.json, N-416, equals "
        "the recalculated threshold; d1_c3.py mode_rulethr gives a one-candidate cell the same arrays and random key "
        "under both rules).",
        "Ledger question (i): one fixed threshold per cell (`judgements.json` `V4.cells.<cell>.threshold`, a number; equal "
        "to `d1_rulethr.json` `cells.<cell>.threshold_stored`), while real coverage varies over the nine k. So both "
        "triangles show coverage − fixed threshold.",
        "Ledger question (ii): not shared within the coil groups (stored values, raw; not printed in the figure). Two "
        f"candidates: fixed thresholds {rng(a[2]['thr'])}, recalculated thresholds {rng(a[2]['thr_r'])}, real-data-rule "
        f"coverage {rng(a[2]['tcr'])}. Four candidates: fixed thresholds {rng(a[4]['thr'])}, recalculated thresholds "
        f"{rng(a[4]['thr_r'])}, real-data-rule coverage {rng(a[4]['tcr'])}. The real-data-rule coverage is shared within an "
        "acquisition type × candidate set (asserted), so each cell keeps its own value.",
        "The per-cell coverages, thresholds and test-data coverages have no per-cell ledger row (as in v10 and v11); in the "
        "figure they are drawn only, never printed. The after-the-verdict values are printed per cell in Supplementary "
        "Table S5 by this script (section below); their extremes per coil group are the ledger rows of Section 4.4.",
        "Tracks: TE and TR are the stage-3 card values (midpoints of the recorded ranges in `decomp/s4/plan.py` ONH; for "
        "NOT2ING-89 and NOT2ING-98 the TR of `decomp/stage1/analyze3.py` NOT2ING_FIXED), exactly as analyze3.py x11 built "
        "the cards; the ranges are not drawn (OXF4GEP TE 68.3–71.0 ms, NOT1ACH TE 70.0–71.2 ms and TR 9.81–10.62 s, "
        "NOT2ING-98 TR 4.031–4.378 s, NOT2ING-89 TR 4.004–4.031 s). Their differences reproduce cards.json for all 28 card "
        "pairs and Table S5 for all its rows (asserted). Scheme letters: plan.py «S» -> A (the three Siemens acquisitions), "
        "«V» -> B (Philips and GE); within each letter the cards record the same scheme.",
        "Receive coil: filled chip = given to the judge (plan.py ONH: OXF1PRI 32, OXF2PRI 64); open chip = recorded in every "
        "sidecar but not given (OXF3TRI 32, N-428; read after the verdict); – = no coil field in any sidecar. The sidecar "
        "records (rv_4a_coil.json) were read after the verdict: fix 1 (verifier finding 7) grades them in the figure's "
        "coil key «(sidecar records, read after the verdict)» and in the caption «open chip and dash, sidecar records "
        "read after the verdict».",
        v["s5_check"] + ".",
        "The NDI cell (OXF1PRI~OXF2PRI, unblinded check) is drawn as two small cells in the footnote line with the same "
        "colours, outline and ✓ rule, and no number: its test-data coverage under the real-data rule has no ledger row "
        "(task: «print numbers only with ledger rows; otherwise Section 4.6»).",
    ]


def table_s5_extra(v):
    L = ["## Supplementary Table S5 after the verdict (printed by this script, not in the figure)", "",
         f"Mode written to F6_supp_edits.json: «{S5_MODE}». Rounded half up to two decimals (the rule of N-382 – N-389). "
         "Each value comes from the ledger copy of `remeasure/explore/d1/d1_rulethr.json` "
         f"({v['D1'].sha[:12]}); «ledger ID» lists a row on the same field where one exists. Exact ties printed half up: "
         f"{len(S5_NEW_TIES)}" + (" (" + "; ".join(S5_NEW_TIES) + ")" if S5_NEW_TIES else "") + ".", "",
         "| element | field | value | printed | ledger ID |", "|---|---|---|---|---|"]
    for t in TABLED:
        L.append(f"| {cell_md(t['element'])} | `{cell_md(t['field'])}` | {cell_md(short(t['value']))} | {t['printed']} | "
                 f"{' · '.join(t['ids']) if t['ids'] else 'no ledger row'} |")
    return L + [""]


if __name__ == "__main__":
    v = read_values()                                   # all data guards; nothing is written before they pass
    print(v["s5_check"])
    print("order by TE:", " · ".join(f"{k} {v['acq'][k]['te']:g}" for k in v["order"]))
    vals5 = s5_values(v)                                # Table S5 after-the-verdict values and their guards (once)
    variants = {m: s5_edits(v, m, vals5) for m in ("columns", "part_e")}
    fig = draw6(v)
    parallel_ok(fig, v["_labels"], CS)
    audit = finish(fig, 6, NAME, CAPTION, v["srcs_judge"], notes6(v), table_s5_extra(v))
    for m, ed in variants.items():
        (OUT / f"F6_supp_edits_{m}.json").write_text(json.dumps(ed, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    (OUT / "F6_supp_edits.json").write_text(json.dumps(variants[S5_MODE], ensure_ascii=False, indent=1) + "\n",
                                            encoding="utf-8")
    print(f"supplement edits: {S5_MODE} ({len(variants[S5_MODE])} edits; old texts unique in {SUPP_MD.name})")
    cvd = FC.cvd_sheet(OUT / f"{NAME}.png", OUT / "cvd" / "F6_cvd.png")
    print("cvd sheet:", cvd.name, hashlib.sha256(cvd.read_bytes()).hexdigest()[:16])
    print("ledger copies read:", len(READS))
    for rel, sha in sorted(READS.items()):
        print(f"  read {rel} {sha[:12]}")
    print(audit)
