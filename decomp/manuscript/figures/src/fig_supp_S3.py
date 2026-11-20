"""Supplementary Fig. S3 — S4 calibration audit, drawn from stored values only (no template mark).

  (a) criterion (i): coverage of the central 90 % intervals by validation stratum and component, C3 and C7,
      against nominal 0.90, the pass threshold 0.90 - 3 sqrt(0.09 / n) and the conservative flag above 0.95;
      grey bands = range over the five ensemble members and three components (criterion iv, reported only)
  (b) criterion (ii): oracle NPE against the MCMC reference posterior — C along u_relax per case
      (TR only and TE and TR), C for ODI and C2ST accuracy per case, each median against its limit
  (c) criterion (v): per-stratum false-alarm rate against 0.02; injection tests (descriptive)
  (d) predictions P1–P6, C3 and C7: filled = met, open = not met, with the values the report prints

Sources, read at run time from the S4 result commit 1f156ed of the read-only repository (git show):
  decomp/results/s4/REPORT.md  §2 lines 11–59 · §3 lines 60–64 · §4 lines 66–81 · §7 lines 124–129
  decomp/results/s4/test_c3.json · test_c7.json  layers[stratum].n / coverage / threshold / fail /
                                                 conservative / coverage_members · n_layers_fail · S4_fail_i · P · S4_fail_P
  decomp/results/s4/oracle.json   rows[].kind / converged / C_ref / C_npe / c2st / rhat_max / ess_bulk_min and summary
  decomp/results/s4/miss.json     layers[stratum].false_alarm · false_alarm_max · false_alarm_ok
  decomp/results/s4/inject.json   [cell].n / detected / n_undetected / coverage90_undetected
  decomp/SPEC_S4.md               lines 120, 122, 127, 129, 135, 136 (levels, pass rule, limits, injections)
No new statistic. The per-case values of (b) are the ones the frozen code reduces to the three medians
(decomp/s4/validate.py:654–659); the script recomputes those medians and stops unless they equal oracle.json.
Every number printed in the figure is compared with REPORT.md at its printed precision, and the script stops
(sys.exit) on any disagreement, on any change of the cited report lines between 1f156ed and HEAD, or if a
working-tree JSON differs from 1f156ed.
Run: CUDA_VISIBLE_DEVICES= MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1 <pdenv>/python -B fig_supp_S3.py
"""
import hashlib
import json
import math
import pathlib
import re
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
import figstyle as fs  # noqa: E402

REPO = pathlib.Path("<repo>")
COMMIT = "1f156ed"                     # S4 result commit cited by the manuscript
S4 = "decomp/results/s4"
SPEC = "decomp/SPEC_S4.md"
NAME = "FigS3_S4_calibration"
COMPONENTS, CCOL, CMRK, CMS = fs.COMPONENTS, fs.CCOL, fs.CMRK, fs.CMS
INK, INK2, MUTED, AXIS, SURF = fs.INK, fs.INK2, fs.MUTED, fs.AXIS, fs.SURF
PT = fs.PT_TEXT
# stratum keys of the JSON (decomp/s4/validate.py layer names), the REPORT row names and the figure names
TYPES = [("match", "일치", "Matched"), ("te", "TE 만", "TE only"), ("tr", "TR 만", "TR only"),
         ("tetr", "TE+TR", "TE and TR"), ("tetr_design", "TE+TR+설계", "TE, TR and scheme")]
LEVELS = ("L0", "L1", "L2", "L3")
INJ = [("bscale", "b-value ±5 %"), ("dpar3.4", "Diffusivity 3.4"), ("neff400", "Neff 400"), ("snr10", "SNR 10")]
DIAMOND, DMS = "D", 3.6                 # neutral pass/fail mark (as Fig. 5): filled = met, open = not met
DOT = "#9d9b94"                         # per-case values in (b)
MEM_C = "#e3e1da"                       # band: range over the five ensemble members in (a)


def stop(msg):
    sys.exit(f"STOP fig S3: {msg}")


def git(*args):
    """Read-only git commands only (show, rev-parse, hash-object without -w)."""
    r = subprocess.run(["git", "-C", str(REPO), *args], capture_output=True)
    if r.returncode:
        stop(f"git {' '.join(args)}: {r.stderr.decode('utf-8', 'replace').strip()}")
    return r.stdout.decode("utf-8")


def show(path, rev=COMMIT):
    return git("show", f"{rev}:{path}")


def f3(v):
    return f"{v:.3f}"


def eq(got, printed, what):
    """A value drawn or printed in the figure must equal the number REPORT.md prints, at its precision."""
    nd = len(printed.split(".")[1]) if "." in printed else 0
    s = f"{got:.{nd}f}"
    if s != printed:
        stop(f"{what}: stored {got!r} -> {s} but REPORT.md prints {printed}")


def match(line_no, pattern, lines, what):
    m = re.fullmatch(pattern, lines[line_no - 1])
    if not m:
        stop(f"REPORT.md:{line_no} ({what}) does not read as expected: {lines[line_no - 1]!r}")
    return m


# ============================================================================ sources
def load():
    rep_txt = show(f"{S4}/REPORT.md")
    head_txt = show(f"{S4}/REPORT.md", "HEAD")
    rep, head = rep_txt.split("\n"), head_txt.split("\n")
    for a, b in ((11, 81), (124, 129)):            # the cited lines must not have moved or changed since 1f156ed
        if rep[a - 1:b] != head[a - 1:b]:
            stop(f"REPORT.md lines {a}–{b} differ between {COMMIT} and HEAD")
    js = {}
    for f in ("test_c3.json", "test_c7.json", "oracle.json", "miss.json", "inject.json"):
        p = f"{S4}/{f}"
        if git("hash-object", str(REPO / p)).strip() != git("rev-parse", f"{COMMIT}:{p}").strip():
            stop(f"working-tree {p} differs from {COMMIT}")
        js[f] = json.loads(show(p))
    spec = show(SPEC).split("\n")
    if spec != show(SPEC, "HEAD").split("\n"):
        stop("SPEC_S4.md differs between the result commit and HEAD")
    return rep, js, spec


# ============================================================================ checks: (i) coverage
def check_coverage(rep, js, spec):
    if not rep[10].startswith("## 2. 기준 (i) 보정"):
        stop(f"REPORT.md:11 is not §2: {rep[10]!r}")
    if "50 · 80 · 90 · 95 %" not in spec[119]:
        stop("SPEC_S4.md:120 levels")
    ci = (0.50, 0.80, 0.90, 0.95)                     # decomp/s4/validate.py:28 CI, as SPEC_S4.md:120
    i90 = ci.index(0.90)
    if "경험 포함률 ≥ 명목 − 3·√(p(1 − p)/n)" not in spec[121] or "과보정(명목 + 0.05 초과)" not in spec[121]:
        stop("SPEC_S4.md:122 pass rule")
    cons = 0.90 + 0.05
    out = {}
    names = {kr: k for k, kr, _ in TYPES}
    for arm, head_no, first in (("c3", 12, 15), ("c7", 36, 39)):
        t = js[f"test_{arm}.json"]
        m = match(head_no, rf"### {arm.upper()} — 실패 층 \*\*(\d+) / (\d+)\*\* · \(i\) 전체 실패\(10 층 이상\) = \*\*(True|False)\*\*",
                  rep, f"{arm} header")
        n_fail = sum(v["fail"] for v in t["layers"].values())
        if not (int(m.group(1)) == n_fail == t["n_layers_fail"] and int(m.group(2)) == len(t["layers"]) == 20
                and (m.group(3) == "True") == t["S4_fail_i"]):
            stop(f"{arm}: failed strata {n_fail} / {len(t['layers'])} vs REPORT {m.groups()}")
        rows = {}
        for ln in range(first, first + 20):
            r = match(ln, r"\| (.+?) · (L[0-3]) \| (\d+) \| (\S+) \| ([0-9.]+) · ([0-9.]+) · ([0-9.]+) \| ([0-9.]+) \| "
                          r"([0-9.]+) – ([0-9.]+) \|", rep, f"{arm} stratum row")
            key = f"{names[r.group(1)]}|{r.group(2)}"
            L = t["layers"][key]
            if int(r.group(3)) != L["n"]:
                stop(f"{arm} {key}: n {L['n']} vs REPORT {r.group(3)}")
            if r.group(4) != "✓" or L["fail"]:
                stop(f"{arm} {key}: fail flag {L['fail']} vs REPORT {r.group(4)}")
            cov = [L["coverage"][c][i90] for c in range(3)]
            for c in range(3):
                eq(cov[c], r.group(5 + c), f"{arm} {key} coverage90 {COMPONENTS[c]}")
                if L["conservative"][c][i90] != (cov[c] > cons):
                    stop(f"{arm} {key}: conservative flag")
            thr = 0.90 - 3 * math.sqrt(0.90 * 0.10 / L["n"])
            if abs(thr - L["threshold"][i90]) > 1e-12:
                stop(f"{arm} {key}: threshold {L['threshold'][i90]} vs rule {thr}")
            eq(thr, r.group(8), f"{arm} {key} threshold")
            # pass at every level and component <=> not failed (decomp/s4/validate.py:309)
            if any(L["coverage"][c][j] < L["threshold"][j] for c in range(3) for j in range(4)):
                stop(f"{arm} {key}: a cell lies below its threshold")
            mem = [L["coverage_members"][k][c][i90] for k in range(len(L["coverage_members"])) for c in range(3)]
            if len(L["coverage_members"]) != 5:
                stop(f"{arm} {key}: {len(L['coverage_members'])} members")
            eq(min(mem), r.group(9), f"{arm} {key} member min")
            eq(max(mem), r.group(10), f"{arm} {key} member max")
            rows[key] = dict(cov=cov, thr=thr, mem=(min(mem), max(mem)), passed=not L["fail"])
        if len(rows) != 20:
            stop(f"{arm}: {len(rows)} strata rows")
        out[arm] = dict(rows=rows, n_fail=n_fail)
    return out, cons


# ============================================================================ checks: (ii) prior leakage
def check_oracle(rep, js, spec):
    o = js["oracle.json"]
    if not rep[59].startswith("## 3. 기준 (ii) 사전 누출"):
        stop(f"REPORT.md:60 is not §3: {rep[59]!r}")
    m1 = match(61, r"- 참조 경우 (\d+) · 수렴 못 해 뺀 경우 \*\*(\d+)\*\*\(([0-9.]+)% · 문턱 (\d+) %\) · .*", rep, "oracle n")
    m2 = match(62, r"- \(a\) TR 만 · TE\+TR 층의 \|C_npe − C_ref\|\(u_relax\) 중앙값 \*\*([0-9.]+)\*\*\(≤ ([0-9.]+)\) · "
                   r"\(b\) C_ODI 중앙값 \*\*([0-9.]+)\*\*\(≥ ([0-9.]+)\) · \(c\) C2ST 정확도 중앙값 \*\*([0-9.]+)\*\*\(≤ ([0-9.]+)\) · "
                   r"NPE 표본 부족 (\d+)", rep, "oracle medians")
    m3 = match(63, r"- \*\*\(ii\) 통과 = (True|False)\*\* · 실패면 S4 전체 실패", rep, "oracle pass")
    m4 = match(64, r"- 수렴 진단: R̂ 최대 ([0-9.]+) · bulk ESS 최소 (\d+)", rep, "oracle diagnostics")
    if "(a) TR 만 · TE+TR 층에서 신탁 NPE 의 C_{u_relax} 와 참조 C 의 차 ≤ 0.10(중앙값) (b) C_ODI ≥ 0.9 (c) C2ST 정확도 중앙값 ≤ 0.65" \
            not in spec[128]:
        stop("SPEC_S4.md:129 criterion (ii)")
    if "획득 유형마다 100 경우" not in spec[126]:
        stop("SPEC_S4.md:127 reference cases")
    lim = dict(a=float(m2.group(2)), b=float(m2.group(4)), c=float(m2.group(6)))
    if lim != dict(a=0.10, b=0.9, c=0.65):
        stop(f"limits {lim}")
    rows = o["rows"]
    n_ex = sum(not r["converged"] for r in rows)
    if not (int(m1.group(1)) == len(rows) == o["n"] and int(m1.group(2)) == n_ex == o["n_excluded_not_converged"]):
        stop("oracle case counts")
    eq(100 * o["frac_excluded"], m1.group(3), "oracle excluded %")
    if int(m1.group(4)) != round(100 * o["mcmc"]["max_excluded"]):
        stop("oracle exclusion limit")
    kinds = {k: sum(r["kind"] == k for r in rows) for k, _, _ in TYPES}
    if set(kinds.values()) != {100}:
        stop(f"oracle cases per type {kinds}")
    # per-case values exactly as the frozen code reduces them (decomp/s4/validate.py:654–659)
    sel = [r for r in rows if "C_npe" in r]
    pa = [(r["C_ref"][3], r["C_npe"][3]) for r in sel if r["kind"] in ("tr", "tetr")
          and np.isfinite(r["C_npe"][3]) and np.isfinite(r["C_ref"][3])]
    if len(pa) != kinds["tr"] + kinds["tetr"]:
        stop(f"criterion (ii)(a) uses {len(pa)} cases, not the {kinds['tr'] + kinds['tetr']} TR-only and TE+TR cases")
    da = [abs(n - c) for c, n in pa]
    codi = [r["C_npe"][2] for r in sel]
    c2 = [r["c2st"] for r in sel]
    n_short = sum(1 for r in rows if r["converged"] and "C_npe" not in r)
    med = dict(a=float(np.median(da)), b=float(np.median(codi)), c=float(np.median(c2)))
    for k, field in (("a", "a_median_abs_diff_C_u_relax"), ("b", "b_median_C_ODI"), ("c", "c_median_c2st")):
        if med[k] != o[field]:
            stop(f"recomputed median {k} {med[k]} != oracle.json {field} {o[field]}")
    eq(med["a"], m2.group(1), "(ii)(a) median")
    eq(med["b"], m2.group(3), "(ii)(b) median")
    eq(med["c"], m2.group(5), "(ii)(c) median")
    if not (int(m2.group(7)) == n_short == o["n_npe_short"] and len(sel) == len(rows) - n_ex - n_short):
        stop("NPE short count")
    passed = med["a"] <= lim["a"] and med["b"] >= lim["b"] and med["c"] <= lim["c"]
    if not (passed == o["pass"] == (m3.group(1) == "True") and o["S4_fail_ii"] == (not passed)):
        stop("(ii) pass")
    eq(max(r["rhat_max"] for r in rows), m4.group(1), "R-hat max")
    eq(min(r["ess_bulk_min"] for r in rows), m4.group(2), "ESS min")
    return dict(pairs=pa, codi=codi, c2st=c2, med=med, lim=lim, n=len(rows), n_ex=n_ex, n_sel=len(sel), passed=passed)


# ============================================================================ checks: (v) misspecification, injections
def check_miss(rep, js, spec):
    mi, inj = js["miss.json"], js["inject.json"]
    if not rep[123].startswith("## 7. 기준 (v) 오지정 점수"):
        stop(f"REPORT.md:124 is not §7: {rep[123]!r}")
    m = match(125, r"- 층별 99 백분위 문턱\(앞 반\) · 뒤 반 오경보율 최대 \*\*([0-9.]+)\*\*\(≤ ([0-9.]+) → (True|False)\) · "
                   r"범위 ([0-9.]+) – ([0-9.]+)", rep, "false alarms")
    if "층별 오경보율 ≤ 2 %" not in spec[134]:
        stop("SPEC_S4.md:135 limit")
    lim = float(m.group(2))
    if lim != 0.02:
        stop(f"false-alarm limit {lim}")
    fa = {k: v["false_alarm"] for k, v in mi["layers"].items()}
    want = {f"{k}|{lv}" for k, _, _ in TYPES for lv in LEVELS}
    if set(fa) != want:
        stop("miss.json strata")
    if max(fa.values()) != mi["false_alarm_max"]:
        stop("false_alarm_max")
    eq(max(fa.values()), m.group(1), "false-alarm max")
    eq(min(fa.values()), m.group(4), "false-alarm min")
    eq(max(fa.values()), m.group(5), "false-alarm range max")
    ok = max(fa.values()) <= lim
    if not (ok == mi["false_alarm_ok"] == (m.group(3) == "True")):
        stop("(v) verdict")
    for s in ("SNR 10", "d∥ 3.4", "b 값 배율 ±5 %", "N_eff 400", "각 1,000 경우"):
        if s not in spec[135]:
            stop(f"SPEC_S4.md:136 injection {s!r}")
    injd = []
    for ln, (key, label) in zip(range(126, 130), INJ):
        r = match(ln, r"- 주입 «(.+?)»: 탐지율 ([0-9.]+) · 탐지 안 된 (\d+) 경우의 90 % 포함률 NDI · FISO · ODI "
                      r"([0-9.]+) · ([0-9.]+) · ([0-9.]+)", rep, "injection")
        if r.group(1) != key:
            stop(f"injection order {r.group(1)} != {key}")
        v = inj[key]
        if v["n"] != 1000 or int(r.group(3)) != v["n_undetected"]:
            stop(f"injection {key} counts")
        eq(v["detected"], r.group(2), f"{key} detected")
        for c in range(3):
            eq(v["coverage90_undetected"][c], r.group(4 + c), f"{key} coverage {COMPONENTS[c]}")
        injd.append(dict(key=key, label=label, det=v["detected"], cov=v["coverage90_undetected"]))
    return dict(fa=fa, lim=lim, ok=ok, inj=injd)


# ============================================================================ checks: P1–P6
def check_p(rep, js):
    if not rep[65].startswith("## 4. 사전 예측 P1 – P6"):
        stop(f"REPORT.md:66 is not §4: {rep[65]!r}")
    B = r"\*\*(True|False)\*\*"
    out = {}
    for arm, h in (("c3", 67), ("c7", 75)):
        P = js[f"test_{arm}.json"]["P"]
        hm = match(h, rf"### {arm.upper()} — P1 · P2 중 하나라도 깨지면 S4 전체 실패 = {B}", rep, f"{arm} P header")
        m1 = match(h + 1, rf"- P1\(일치 · L0 · 세 성분 C ≥ 0\.75\): NDI ([0-9.]+) · FISO ([0-9.]+) · ODI ([0-9.]+) → {B}", rep, "P1")
        m2 = match(h + 2, rf"- P2\(모든 획득 유형 · L0 – L2 · C_ODI ≥ 0\.75\): 15 층의 최솟값 ([0-9.]+) → {B}", rep, "P2")
        m3 = match(h + 3, rf"- P3\(TR 만 · L0 의 C_u_relax ≤ 일치 · L0 의 C_u_ref − 0\.2 · C_ODI 차 ≤ 0\.1\): ([0-9.]+) 대 ([0-9.]+) · "
                          rf"ODI ([0-9.]+) 대 ([0-9.]+) → {B}", rep, "P3")
        m4 = match(h + 4, rf"- P4\(TE 만 · L0 의 C_NDI < 일치 · L0\): ([0-9.]+) 대 ([0-9.]+) → {B}", rep, "P4")
        m5 = match(h + 5, rf"- P5\(일치 · L1 의 σ_NDI 증가가 0\.397 × s_rms 의 ± 50 %\): σ L1\(ρ 0\) ([0-9.]+) · L0 ([0-9.]+) · "
                          rf"Δσ ([0-9.]+) · 기대 ([0-9.]+) → {B}", rep, "P5")
        m6 = match(h + 6, rf"- P6\(L3 에서 u_relax 의 NDI 성분 부호 반전 비율 > 0\.5 · 대수 · k 무관\): ([0-9.]+) → {B}", rep, "P6")
        p1 = P["P1"]["median_of_k"]
        for i, c in enumerate(COMPONENTS):
            eq(p1[c], m1.group(1 + i), f"{arm} P1 {c}")
        p2 = P["P2"]["median_of_k"]
        if len(p2) != 15:
            stop(f"{arm} P2 strata {len(p2)}")
        eq(min(p2.values()), m2.group(1), f"{arm} P2 lowest")
        p3 = P["P3"]["median_of_k"]
        for i, k in enumerate(("C_u_relax_TR", "C_u_ref_match", "C_ODI_TR", "C_ODI_match")):
            eq(p3[k], m3.group(1 + i), f"{arm} P3 {k}")
        p4 = P["P4"]["median_of_k"]
        eq(p4["C_NDI_TE"], m4.group(1), f"{arm} P4 TE")
        eq(p4["C_NDI_match"], m4.group(2), f"{arm} P4 matched")
        p5 = P["P5"]
        eq(p5["median_of_k"]["sd_L1_rho0"], m5.group(1), f"{arm} P5 sd L1")
        eq(p5["median_of_k"]["sd_L0"], m5.group(2), f"{arm} P5 sd L0")
        eq(p5["dsigma"], m5.group(3), f"{arm} P5 dsigma")
        eq(p5["expected"], m5.group(4), f"{arm} P5 expected")
        eq(P["P6"]["frac_sign_flip"], m6.group(1), f"{arm} P6")
        holds = {f"P{i}": P[f"P{i}"]["holds"] for i in range(1, 7)}
        for i, mm in enumerate((m1, m2, m3, m4, m5, m6), 1):
            if holds[f"P{i}"] != (mm.groups()[-1] == "True"):
                stop(f"{arm} P{i}: stored {holds[f'P{i}']} vs REPORT {mm.groups()[-1]}")
        fail_p = not (holds["P1"] and holds["P2"])
        if not (fail_p == js[f"test_{arm}.json"]["S4_fail_P"] == (hm.group(1) == "True")):
            stop(f"{arm} P1/P2 overall")
        vals = {"P1": " · ".join(f3(p1[c]) for c in COMPONENTS),
                "P2": f"lowest {f3(min(p2.values()))}",
                "P3": f"{f3(p3['C_u_relax_TR'])} vs {f3(p3['C_u_ref_match'])}",
                "P4": f"{f3(p4['C_NDI_TE'])} vs {f3(p4['C_NDI_match'])}",
                "P5": f"Δσ {f3(p5['dsigma'])} vs {f3(p5['expected'])}",
                "P6": f3(P["P6"]["frac_sign_flip"])}
        out[arm] = dict(holds=holds, vals=vals)
    return out


# ============================================================================ drawing
def dmark(ax, x, y, met, ms=DMS, color=INK, mew=0.9, **kw):
    ax.plot([x], [y], ls="none", marker=DIAMOND, ms=ms, mfc=color if met else SURF, mec=color,
            mew=mew, zorder=kw.pop("zorder", 5), clip_on=kw.pop("clip_on", False), **kw)


def cmark(ax, x, y, c, filled=True, zorder=4, ms=None):
    ax.plot([x], [y], ls="none", marker=CMRK[c], ms=ms or CMS[c], mfc=CCOL[c] if filled else SURF,
            mec=CCOL[c] if not filled else SURF, mew=0.9 if not filled else 0.4, zorder=zorder, clip_on=False)


def jitter(n, width):
    """Deterministic spread for a strip (golden-ratio sequence); positions only, values untouched."""
    return (((np.arange(n) * 0.6180339887498949) % 1.0) - 0.5) * width


def draw(cov, cons, orc, mis, pp):
    W, H = 190.0, 150.0
    fig = fs.new_fig(W, H)
    cv = fs.canvas(fig)

    # ------------------------------------------------------------------ (a) coverage by stratum
    XA0, XA1 = 23.0, 160.0
    RH, YC7 = 21.0, 93.0
    YC3 = YC7 + RH + 3.0
    xs, groups = {}, []
    x = 0.0
    for k, _, _ in TYPES:
        c0 = x
        for lv in LEVELS:
            xs[f"{k}|{lv}"] = x
            x += 1.0
        groups.append((c0, x - 1.0))
        x += 0.75
    xmax = x - 0.75 - 1.0
    off = {"NDI": -0.25, "FISO": 0.0, "ODI": 0.25}
    YLO, YHI = 0.80, 0.965
    axs = {}
    for arm, y0 in (("c3", YC3), ("c7", YC7)):
        a = fs.ax_mm(fig, XA0, y0, XA1 - XA0, RH)
        axs[arm] = a
        a.set_xlim(-0.6, xmax + 0.6)
        a.set_ylim(YLO, YHI)
        rows = cov[arm]["rows"]
        # fail zone below each stratum's own threshold (0.90 - 3 sqrt(0.09 / n), n = cases kept in that stratum)
        for key, xx in xs.items():
            thr = rows[key]["thr"]
            a.fill_between([xx - 0.5, xx + 0.5], YLO, thr, color=fs.ZONE_PART, lw=0, zorder=0)
            a.plot([xx - 0.5, xx + 0.5], [thr, thr], color=INK2, lw=0.6, ls=(0, (2.2, 1.4)), zorder=1)
        a.axhline(0.90, color=AXIS, lw=0.6, zorder=1.0)
        a.axhline(cons, color=AXIS, lw=0.5, ls=(0, (1, 1.6)), zorder=1.0)
        for c0, _ in groups[1:]:
            a.axvline(c0 - 0.875, color=fs.GRID, lw=0.5, zorder=0.5)
        for key, xx in xs.items():
            r = rows[key]
            a.add_patch(Rectangle((xx - 0.40, r["mem"][0]), 0.80, r["mem"][1] - r["mem"][0], facecolor=MEM_C,
                                  edgecolor="none", lw=0, zorder=0.9))
            for i, c in enumerate(COMPONENTS):
                cmark(a, xx + off[c], r["cov"][i], c, filled=r["passed"])
        a.set_yticks([0.80, 0.85, 0.90, 0.95])
        a.set_yticklabels(["0.80", "0.85", "0.90", "0.95"])
        a.tick_params(axis="x", length=0)
        if arm == "c7":
            a.set_xticks(list(xs.values()))
            a.set_xticklabels([k.split("|")[1] for k in xs])
            a.tick_params(axis="x", pad=1.5)
        else:
            a.set_xticks([])
        a.spines["bottom"].set_visible(arm == "c7")
        cv.text(2.0, y0 + RH - 2.0, arm.upper(), fontsize=PT, fontweight="bold", color=INK, ha="left", va="center")
        cv.text(2.0, y0 + RH - 5.0, f"{cov[arm]['n_fail']} of {len(rows)}\nstrata\nfailed", fontsize=PT, color=INK2,
                ha="left", va="top", linespacing=1.1)
        thr_txt = {f3(rows[k]["thr"]) for k in rows}          # every stratum prints 0.871 (REPORT.md table)
        if len(thr_txt) != 1:
            stop(f"{arm}: thresholds print as {thr_txt}")
        tr = a.get_yaxis_transform()
        for yv, txt, va, col in ((cons, f"{cons:.2f} conservative", "center", INK2),
                                 (0.90, "0.90 nominal", "center", INK2),
                                 (0.8745, f"{thr_txt.pop()} pass threshold", "top", INK2),
                                 (0.825, "fail zone", "center", MUTED)):
            a.text(1.01, yv, txt, transform=tr, fontsize=PT, color=col, va=va, ha="left", clip_on=False)
    a7 = axs["c7"]
    for (c0, c1), (_, _, lab) in zip(groups, TYPES):
        a7.annotate(lab, ((c0 + c1) / 2, 0), xycoords=("data", "axes fraction"), xytext=(0, -12.0),
                    textcoords="offset points", ha="center", va="top", fontsize=PT, color=INK, annotation_clip=False)
    cv.text(13.0, YC7 + RH + 1.5, "Coverage of the central 90 % interval", rotation=90, fontsize=PT, color=INK,
            ha="center", va="center")
    LET_A = YC3 + RH + 4.5
    fs.panel(axs["c3"], "a", "Calibration, criterion (i), true k", dx_mm=4.0 - XA0, dy_mm=LET_A - (YC3 + RH))
    ky = LET_A + 1.6
    kx = 80.0
    for c, dx in zip(COMPONENTS, (0.0, 9.6, 21.0)):
        cmark(cv, kx + dx, ky, c, filled=True)
        cv.text(kx + dx + 1.7, ky, c, fontsize=PT, color=INK, ha="left", va="center")
    cv.text(kx + 29.0, ky, "ensemble; filled = stratum passed", fontsize=PT, color=INK, va="center")
    cv.add_patch(Rectangle((kx + 71.0, ky - 1.6), 2.4, 3.2, facecolor=MEM_C, edgecolor="none", lw=0))
    cv.text(kx + 74.8, ky, "range of the five members", fontsize=PT, color=INK, va="center")

    # ------------------------------------------------------------------ (b) prior leakage
    LET_B = 78.5                         # common baseline of the letters (b) and (c)
    sb = fs.ax_mm(fig, 19.0, 42.0, 27.0, 27.0)
    pa = np.array(orc["pairs"])
    sb.plot([0.2, 1.0], [0.2, 1.0], color=AXIS, lw=0.6, zorder=1)
    sb.plot(pa[:, 0], pa[:, 1], ls="none", marker="o", ms=1.7, mfc=DOT, mec="none", alpha=0.75, zorder=2)
    sb.set_xlim(0.2, 1.0)
    sb.set_ylim(0.2, 1.0)
    for setter in (sb.set_xticks, sb.set_yticks):
        setter([0.2, 0.6, 1.0])
    sb.set_xticklabels(["0.2", "0.6", "1"])
    sb.set_yticklabels(["0.2", "0.6", "1"])
    sb.set_xlabel("C, reference (MCMC)", labelpad=1.5)
    sb.set_ylabel("C, oracle NPE", labelpad=1.5)
    sb.text(0.04, 0.97, f"u_relax, n = {len(pa)}", transform=sb.transAxes, fontsize=PT, color=INK, ha="left",
            va="top")
    sb.text(0.40, 0.235, f"median |ΔC| {f3(orc['med']['a'])}: {'met' if orc['med']['a'] <= orc['lim']['a'] else 'not met'}",
            transform=sb.transAxes, fontsize=PT, color=INK, ha="left", va="center")
    sb.text(0.40, 0.09, f"limit ≤ {orc['lim']['a']:.2f}", transform=sb.transAxes, fontsize=PT, color=INK2,
            ha="left", va="center")
    fs.panel(sb, "b", "Prior leakage, criterion (ii)", dx_mm=4.0 - 19.0, dy_mm=LET_B - 69.0)

    def strip(y0, vals, med, lim, side, xlim, ticks, title, fmt):
        s = fs.ax_mm(fig, 58.0, y0, 43.0, 5.5)
        v = np.asarray(vals)
        s.plot(v, 0.5 + jitter(len(v), 0.66), ls="none", marker="o", ms=1.5, mfc=DOT, mec="none", alpha=0.7,
               zorder=2)
        s.axvline(lim, color=INK2, lw=0.6, ls=(0, (2.2, 1.4)), zorder=1)
        dmark(s, med, 0.5, med >= lim if side == "ge" else med <= lim, zorder=6, clip_on=True)
        s.set_xlim(*xlim)
        s.set_ylim(0, 1)
        s.set_yticks([])
        s.spines["left"].set_visible(False)
        s.set_xticks(ticks)
        s.set_xticklabels([fmt(t) for t in ticks])
        s.text(0.0, 1.0 + 3.6 / 5.5, title, transform=s.transAxes, fontsize=PT, color=INK, ha="left", va="bottom")
        s.text(0.0, 1.0 + 0.8 / 5.5, f"median {f3(med)}, limit {'≥' if side == 'ge' else '≤'} {lim:g}",
               transform=s.transAxes, fontsize=PT, color=INK2, ha="left", va="bottom")
        return s

    strip(64.0, orc["codi"], orc["med"]["b"], orc["lim"]["b"], "ge", (0.88, 1.0), [0.88, 0.92, 0.96, 1.0],
          f"C for ODI, oracle NPE, n = {orc['n_sel']}", lambda t: f"{t:.2f}")
    strip(44.5, orc["c2st"], orc["med"]["c"], orc["lim"]["c"], "le", (0.5, 0.9), [0.5, 0.6, 0.7, 0.8, 0.9],
          f"C2ST accuracy, n = {orc['n_sel']}", lambda t: f"{t:.1f}")

    # ------------------------------------------------------------------ (c) misspecification and injections
    fa = mis["fa"]
    fb = fs.ax_mm(fig, 121.0, 66.5, 66.0, 10.0)
    stacks = {}
    for key in [f"{k}|{lv}" for k, _, _ in TYPES for lv in LEVELS]:
        stacks.setdefault(round(fa[key], 6), []).append(key)
    for v, keys in stacks.items():
        for j, _ in enumerate(keys):
            dmark(fb, v, j, v <= mis["lim"], ms=2.5, mew=0.7, zorder=4)
    fb.axvline(mis["lim"], color=INK2, lw=0.6, ls=(0, (2.2, 1.4)), zorder=1)
    fb.set_xlim(0.0, 0.04)
    fb.set_ylim(-0.6, 4.9)
    fb.set_yticks([])
    fb.spines["left"].set_visible(False)
    fb.set_xticks([0, 0.01, 0.02, 0.03, 0.04])
    fb.set_xticklabels(["0", "0.01", "0.02", "0.03", "0.04"])
    fb.set_xlabel("False-alarm rate per stratum (20 strata)", labelpad=1.5)
    over = sorted(k for k, v in fa.items() if v > mis["lim"])
    names = {k: n for k, _, n in TYPES}
    kinds = sorted({k.split("|")[0] for k in over})
    if len(kinds) != 1:
        stop(f"strata above the limit span {kinds}; the label below names one type")
    lvls = [k.split("|")[1] for k in over]
    vmax = max(fa.values())
    if len(stacks[round(vmax, 6)]) != len(over):
        stop("open marks and the strata above the limit differ")
    fb.text(vmax + 0.0012, -0.3, f"{names[kinds[0]]},\n{' and '.join(lvls)}: {f3(vmax)}", fontsize=PT, color=INK,
            ha="left", va="bottom", linespacing=1.1)
    fb.text(mis["lim"] + 0.0007, 4.9, f"limit {mis['lim']:g}: {'met' if mis['ok'] else 'not met'}", fontsize=PT,
            color=INK, ha="left", va="top")
    fs.panel(fb, "c", "Misspecification, criterion (v), C3 arm", dx_mm=108.0 - 121.0, dy_mm=LET_B - 76.5)

    ib = fs.ax_mm(fig, 144.0, 32.0, 43.0, 22.5)
    nrow = len(mis["inj"])
    for i, d in enumerate(mis["inj"]):
        yy = nrow - 1 - i
        for j, c in enumerate(COMPONENTS):
            ib.plot([d["cov"][j]], [yy + (1 - j) * 0.32], ls="none", marker="|", ms=5.0, mew=1.2, color=fs.CCOL[c],
                    zorder=3)                                  # descriptive: no fill (no pass criterion)
    ib.axvline(0.90, color=AXIS, lw=0.6, zorder=0.8)
    ib.set_ylim(-0.6, nrow - 0.4)
    ib.set_xlim(0.55, 0.95)
    ib.set_xticks([0.6, 0.7, 0.8, 0.9])
    ib.set_xticklabels(["0.6", "0.7", "0.8", "0.9"])
    ib.set_yticks(range(nrow))
    ib.set_yticklabels([f"{d['label']}\ndetected {f3(d['det'])}" for d in mis["inj"]][::-1], linespacing=1.05)
    ib.tick_params(axis="y", length=0, pad=2)
    fs.hgrid(ib, axis="x")
    ib.text(-23.0 / 43.0, 1.0 + 1.0 / 22.5, "Injections (descriptive)", transform=ib.transAxes, fontsize=PT,
            color=INK, ha="left", va="bottom")
    ib.set_xlabel("90 % coverage of undetected cases", labelpad=1.5)

    # ------------------------------------------------------------------ (d) predictions P1–P6
    labels = {
        "P1": "Matched, L0: C ≥ 0.75 for NDI, FISO and ODI",
        "P2": "Every acquisition type, L0–L2: C ≥ 0.75 for ODI",
        "P3": "TR only, L0: C along u_relax ≥ 0.2 below matched, ODI within 0.1",
        "P4": "TE only, L0: C for NDI below matched",
        "P5": "Matched: added NDI spread √(σ²L1 − σ²L0) within ±50 % of 0.397·rms(s_T2)",
        "P6": "L3: NDI part of u_relax reverses sign (share > 0.5)",
    }
    db = fs.ax_mm(fig, 4.0, 2.5, 182.0, 20.5)
    db.axis("off")
    db.set_xlim(0, 182)
    db.set_ylim(0, 20.5)
    XL, XC3, XC7 = 9.0, 113.0, 149.0
    ytop = 18.6
    db.text(XL, ytop, "Prediction (P1 and P2 decisive)", fontsize=PT, color=INK2, va="center")
    for arm, xc in (("c3", XC3), ("c7", XC7)):
        db.text(xc, ytop, arm.upper(), fontsize=PT, color=INK, fontweight="bold", va="center", ha="center")
        db.text(xc + 3.0, ytop, "values", fontsize=PT, color=INK2, va="center")
    db.plot([0, 182], [ytop - 1.8, ytop - 1.8], color=AXIS, lw=0.5)
    for i, p in enumerate(("P1", "P2", "P3", "P4", "P5", "P6")):
        y = ytop - 4.3 - i * 2.85
        db.text(0.5, y, p, fontsize=PT, color=INK, fontweight="bold" if p in ("P1", "P2") else "normal", va="center")
        db.text(XL, y, labels[p], fontsize=PT, color=INK, va="center")
        for arm, xc in (("c3", XC3), ("c7", XC7)):
            dmark(db, xc, y, pp[arm]["holds"][p])
            db.text(xc + 3.0, y, pp[arm]["vals"][p], fontsize=PT, color=INK, va="center")
    fs.panel(db, "d", "Predictions fixed before validation", dx_mm=0.0, dy_mm=1.2)
    # one key for every diamond in (b)–(d), in the header row of (d)
    dmark(db, 76.0, ytop, True)
    db.text(77.6, ytop, "met", fontsize=PT, color=INK2, va="center")
    dmark(db, 85.0, ytop, False)
    db.text(86.6, ytop, "not met", fontsize=PT, color=INK2, va="center")
    return fig


def main():
    rep, js, spec = load()
    cov, cons = check_coverage(rep, js, spec)
    orc = check_oracle(rep, js, spec)
    mis = check_miss(rep, js, spec)
    pp = check_p(rep, js)
    print(f"checks ok: (i) failed strata C3 {cov['c3']['n_fail']} C7 {cov['c7']['n_fail']} · (ii) medians "
          f"{f3(orc['med']['a'])} {f3(orc['med']['b'])} {f3(orc['med']['c'])} pass {orc['passed']} · (v) max "
          f"{f3(max(mis['fa'].values()))} ok {mis['ok']} · P C3 {pp['c3']['holds']} · P C7 {pp['c7']['holds']}")
    fig = draw(cov, cons, orc, mis, pp)
    fs.save(fig, NAME)
    for ext in ("pdf", "png"):
        f = HERE.parent / f"{NAME}.{ext}"
        print(f"sha256 {f.name} {hashlib.sha256(f.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
