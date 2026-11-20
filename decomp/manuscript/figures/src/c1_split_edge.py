"""E4 (exploratory): categories 2 and 3 of the (c-1) count, each split into the part at the outer edge of the template brain
mask and the part in or next to its unfilled ventricular gaps (drafts/x15/v12/V12C_S44_RULES.md §1, recorded 2026-11-07
17:13:29 KST, sha256 4b738608…; RULES_RECORD.jsonl). Exploratory: the totals had already been seen in a check outside the
rules (v12b Fig. 7 work, 2026-11-07 03:28 KST); this script fixes that check's definitions and reports them. Verdicts
unchanged (#34); no new effective L; no participant file and no raw image is read.

Inputs: drafts/x15/v12/out/c1_mask_location.json and c1_projections.npz (the per-session template voxels of the warped
masks; npz sha256 = the JSON); the two templates MNI152_T1_1mm_brain_mask and HarvardOxford-sub-maxprob-thr25-1mm, read
again after guard.require() and logged in drafts/x14/EXTERNAL_READ_LOG.jsonl with the functions of c1_mask_location.py
(loaded as a module); their sha256 must equal those recorded in c1_mask_location.json.
Definitions (§1): BM = brain mask != 0; F = binary_fill_holes(BM); G = F and not BM (the unfilled gaps); LV = label 3 or 14.
  category 2 (v12b: not BM, not LV) = 2a (outside F) + 2b (in G);
  category 3 (v12b: BM, distance to the nearest voxel outside BM <= 5 mm, not LV) = 3a (distance to the nearest voxel
  outside F <= 5 mm: the outer edge) + 3b (the rest: within 5 mm only of a gap); distances distance_transform_edt, 1 mm;
  categories 1 and 4 as v12b. Shares per session over its template voxels; S_out = 2a + 3a, S_gap = 2b + 3b.
Guards: categories 1 - 4 recomputed from the re-read templates equal the v12b counts of every session; 2a + 2b = category 2
  and 3a + 3b = category 3 in every session; the rules file and its record; the inputs' sha256.
Report: per scanner and over all sessions, median, min and max of each share (raw and half up to two decimals); over the
  summed voxels of all sessions, the gap-side share of category 2 (2b / 2), of category 3 (3b / 3) and of both
  ((2b + 3b) / (2 + 3)); the name rule of §1 (mostly / partly, per category, from its summed gap-side share).
Writes drafts/x15/v12/out/c1_split.json only.
Run (CPU only):
  PYTHONDONTWRITEBYTECODE=1 CUDA_VISIBLE_DEVICES= DMRI_REPO=<repo> C1_TMP=<scratch>/tmp nice -n 10 \
    <repo>/.venv/bin/python -I -B c1_split_edge.py
"""
import sys

sys.dont_write_bytecode = True

import datetime  # noqa: E402
import importlib.util  # noqa: E402
import io  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import pathlib  # noqa: E402
import statistics  # noqa: E402

_spec = importlib.util.spec_from_file_location("c1_mask_location", pathlib.Path(__file__).resolve().with_name(
    "c1_mask_location.py"))
C = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C)               # functions only (main() is not called); env checks of that module apply
np, ndi = C.np, C.ndi

RULES_C = C.V12 / "V12C_S44_RULES.md"
RULES_C_SHA = "4b73860865207b2235d447d209bbdeabe743def383f9ed7e261dedaa624a668a"
C.RUN_ID = "C1S-" + datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
C.WHAT = "v12c 탐색 · E4 범주 2 · 3 갈래(바깥 가장자리 · 빈 곳) — exploratory · 판정 불변 #34"
SCRIPT = pathlib.Path(__file__).resolve()
OUT_JSON = C.OUTDIR / "c1_split.json"
SIE = C.SIE
PARTS = ("2a", "2b", "3a", "3b")
NAMES = {"2": ("Outside the brain mask", "in its unfilled ventricular gaps"),
         "3": ("Within 5 mm of its edge", "next to these gaps")}


def stop(msg):
    sys.exit(f"STOP C1S: {msg}")


def rules_check():
    if C.sha_bytes(RULES_C.read_bytes()) != RULES_C_SHA:
        stop("V12C_S44_RULES.md differs from the recorded sha256")
    recs = [json.loads(s) for s in C.RULES_REC.read_text(encoding="utf-8").splitlines() if s.strip()]
    hit = [r for r in recs if r.get("sha256") == RULES_C_SHA and str(r.get("파일", "")).endswith("V12C_S44_RULES.md")]
    if len(hit) != 1:
        stop("RULES_RECORD.jsonl does not hold V12C_S44_RULES.md once")
    return hit[0]["시각"]


def stat(vals):
    return dict(median=float(statistics.median(vals)), min=float(min(vals)), max=float(max(vals)), n=len(vals))


def text_of(d):
    return {k: C.half_up(d[k]) for k in ("median", "min", "max")}


def main():
    t0 = C.now()
    rules_t = rules_check()
    if OUT_JSON.exists() and os.environ.get("C1S_OVERWRITE") != "1":
        stop(f"{OUT_JSON} exists (set C1S_OVERWRITE=1 to replace)")
    jp, npz_p = C.OUTDIR / "c1_mask_location.json", C.OUTDIR / "c1_projections.npz"
    jraw = jp.read_bytes()
    J = json.loads(jraw)
    nraw = npz_p.read_bytes()
    if C.sha_bytes(nraw) != J["projections"]["sha256"]:
        stop("c1_projections.npz differs from the sha256 in c1_mask_location.json")
    if J["rules"]["sha256"] != C.RULES_SHA:
        stop("c1_mask_location.json was not made under V12B_S43_RULES.md")
    z = np.load(io.BytesIO(nraw), allow_pickle=False)
    basis = C.G.require()
    bm_raw = C.ext_read(C.T_BM, "MNI152_T1_1mm_brain_mask — E4 갈래(채운 뇌 마스크 F · 빈 곳 G · 바깥 가장자리 거리)", basis)
    ho_raw = C.ext_read(C.T_HO, "HarvardOxford-sub-maxprob-thr25-1mm — E4 갈래(측뇌실 라벨 LV 제외 · 범주 1 대조)", basis)
    tsha = J["grid"]["templates"]["sha256"]
    if C.sha_bytes(bm_raw) != tsha["BM"] or C.sha_bytes(ho_raw) != tsha["HO"]:
        stop("a template differs from the one c1_mask_location.py read")
    bm, ho = C.nii(bm_raw), C.nii(ho_raw)
    if [int(x) for x in bm.shape] != J["grid"]["templates"]["shapes"]["BM"] or \
            float(np.abs(bm.affine - ho.affine).max()) > C.TOL:
        stop("template grid")
    BM = np.asanyarray(bm.dataobj) != 0
    LV = np.isin(np.asanyarray(ho.dataobj), C.LV_VALUES)
    F = ndi.binary_fill_holes(BM)
    Gp = F & ~BM
    d_bm = ndi.distance_transform_edt(BM, sampling=(1.0, 1.0, 1.0))
    d_f = ndi.distance_transform_edt(F, sampling=(1.0, 1.0, 1.0))
    cat = np.full(BM.shape, 4, np.uint8)            # v12b categories, recomputed (guard)
    cat[BM & (d_bm <= C.BAND_MM)] = 3
    cat[~BM] = 2
    cat[LV] = 1
    part = np.zeros(BM.shape, np.uint8)             # 0: categories 1 and 4; 1 2a, 2 2b, 3 3a, 4 3b
    part[(cat == 2) & ~F] = 1
    part[(cat == 2) & Gp] = 2
    part[(cat == 3) & (d_f <= C.BAND_MM)] = 3
    part[(cat == 3) & ~(d_f <= C.BAND_MM)] = 4
    if int(((cat == 2) & (part == 0)).sum()) or int(((cat == 3) & (part == 0)).sum()):
        stop("a category-2 or -3 template voxel is in no part")
    tmpl = dict(n_brain=int(BM.sum()), n_filled=int(F.sum()), n_gap=int(Gp.sum()), n_lv=int(LV.sum()),
                n_lv_in_gap=int((LV & Gp).sum()), n_lv_outside_F=int((LV & ~F).sum()),
                n_template={p: int((part == i + 1).sum()) for i, p in enumerate(PARTS)},
                v12b_n_template=J["template_categories"]["n_template"])
    for c_, n_ in J["template_categories"]["n_template"].items():
        if int((cat == int(c_)).sum()) != n_:
            stop(f"template category {c_}: {int((cat == int(c_)).sum())} voxels, v12b {n_}")
    rows, guards = [], []
    for r in J["sessions"]:
        W = z[f"warped__{r['session']}"].astype(np.int64)
        t = tuple(W.T)
        c = np.bincount(cat[t], minlength=5)
        p = np.bincount(part[t], minlength=5)
        n = len(W)
        if n != r["n_template_voxels"]:
            stop(f"{r['session']}: {n} stored voxels, JSON {r['n_template_voxels']}")
        same = {k: int(c[int(k)]) == r["n_by_category"][k] for k in ("1", "2", "3", "4")}
        g2 = int(p[1] + p[2]) == r["n_by_category"]["2"]
        g3 = int(p[3] + p[4]) == r["n_by_category"]["3"]
        guards.append(dict(session=r["session"], cat_1_4_equal_v12b=all(same.values()), a2_plus_b2_eq_cat2=g2,
                           a3_plus_b3_eq_cat3=g3))
        if not (all(same.values()) and g2 and g3):
            stop(f"{r['session']}: guard failed ({same}, 2a+2b {g2}, 3a+3b {g3})")
        cnt = {pp: int(p[i + 1]) for i, pp in enumerate(PARTS)}
        sh = {pp: cnt[pp] / n for pp in PARTS}
        rows.append(dict(session=r["session"], scanner=r["scanner"], n_template_voxels=n, n_by_part=cnt,
                         n_cat2=r["n_by_category"]["2"], n_cat3=r["n_by_category"]["3"], share_by_part=sh,
                         S_out=sh["2a"] + sh["3a"], S_gap=sh["2b"] + sh["3b"], S_v12b=r["S_share"]))
        if abs((sh["2a"] + sh["2b"] + sh["3a"] + sh["3b"]) - r["S_share"]) > 1e-12:
            stop(f"{r['session']}: S_out + S_gap differs from the v12b S share")
    keys = PARTS + ("S_out", "S_gap")

    def val(r, k):
        return r["share_by_part"][k] if k in PARTS else r[k]
    per_scanner = {dev: {k: stat([val(r, k) for r in rows if r["scanner"] == dev]) for k in keys} for dev in SIE}
    pooled = {k: stat([val(r, k) for r in rows]) for k in keys}
    tot = {pp: sum(r["n_by_part"][pp] for r in rows) for pp in PARTS}
    tot2, tot3 = tot["2a"] + tot["2b"], tot["3a"] + tot["3b"]
    gap = dict(cat2=tot["2b"] / tot2, cat3=tot["3b"] / tot3, both=(tot["2b"] + tot["3b"]) / (tot2 + tot3))
    names = {}
    for c_, k in (("2", "cat2"), ("3", "cat3")):
        word = "mostly" if gap[k] > 0.5 else "partly"
        names[c_] = dict(gap_share=gap[k], rule="summed gap-side share over the 50 sessions > 0.5 -> mostly, else partly",
                         word=word, name=f"{NAMES[c_][0]} ({word} {NAMES[c_][1]})")
    names["both_same_word_under_combined_share"] = ("mostly" if gap["both"] > 0.5 else "partly")
    res = dict(
        what="E4 (exploratory): categories 2 and 3 of the (c-1) count split into the outer edge of the template brain mask "
             "and its unfilled ventricular gaps; totals already seen in a check outside the rules (2026-11-07 03:28 KST); "
             "verdicts unchanged (#34); no new effective L; no participant file read",
        grade="exploratory [탐색]",
        rules=dict(file="drafts/x15/v12/V12C_S44_RULES.md", sha256=RULES_C_SHA, recorded=rules_t, sections="§1"),
        run=dict(id=C.RUN_ID, started=t0, finished=C.now(), script=C.rel(SCRIPT),
                 script_sha256=C.sha_bytes(SCRIPT.read_bytes()), module=C.rel(C.SCRIPT),
                 module_sha256=C.sha_bytes(C.SCRIPT.read_bytes()), guard_basis=basis),
        inputs=dict(c1_mask_location_json=dict(sha256=C.sha_bytes(jraw), run=J["run"]["id"]),
                    c1_projections_npz=dict(sha256=C.sha_bytes(nraw)),
                    templates_sha256=dict(BM=C.sha_bytes(bm_raw), HO=C.sha_bytes(ho_raw)),
                    templates_equal_v12b_reads=True),
        definitions=dict(BM="MNI152_T1_1mm_brain_mask != 0", F="scipy.ndimage.binary_fill_holes(BM)", G="F and not BM",
                         LV="HarvardOxford-sub-maxprob-thr25-1mm value 3 or 14",
                         **{"2a": "category 2 outside F", "2b": "category 2 in G",
                            "3a": "category 3 with distance_transform_edt(F) <= 5 mm (outer edge)",
                            "3b": "category 3, the rest (within 5 mm only of a gap)"},
                         S_out="2a + 3a", S_gap="2b + 3b", denominator="the session's template voxels of value 1",
                         rounding="text: half up to two decimals (Decimal of repr)"),
        template=tmpl,
        sessions=rows, guards=dict(all_pass=True, per_session=guards),
        per_scanner=per_scanner, pooled=pooled,
        totals=dict(n_by_part=tot, n_cat2=tot2, n_cat3=tot3, gap_side_share=gap),
        names=names,
        table_note="the reading \"mostly in the lateral ventricles\" (case ㄴ of the rule) does not change",
        text=dict(per_scanner={dev: {k: text_of(v) for k, v in d.items()} for dev, d in per_scanner.items()},
                  pooled={k: text_of(v) for k, v in pooled.items()},
                  gap_side_share={k: C.half_up(v) for k, v in gap.items()}),
        reads=C.READS)
    OUT_JSON.write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(f"written {OUT_JSON} (sha256 {C.sha_bytes(OUT_JSON.read_bytes())})")
    for dev in SIE:
        print(dev, " · ".join(f"{k} {v['median']:.3f} [{v['min']:.3f}–{v['max']:.3f}]" for k, v in per_scanner[dev].items()))
    print("all", " · ".join(f"{k} {v['median']:.3f} [{v['min']:.3f}–{v['max']:.3f}]" for k, v in pooled.items()))
    print("totals", tot, "gap share", gap, "names", {k: v["name"] for k, v in names.items() if isinstance(v, dict)})
    print("reads", [(r["log"], r["line"]) for r in C.READS])


if __name__ == "__main__":
    main()
