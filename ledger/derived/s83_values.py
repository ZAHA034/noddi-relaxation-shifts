"""Further values from the design-simulation results (no new computation: per-cell metrics copied and counted)."""
import glob, hashlib, json, pathlib, sys

HERE = pathlib.Path(__file__).resolve().parent
OUT = pathlib.Path(sys.argv[1])
cells = {}
for p in sorted(glob.glob(str(OUT / "cells" / "cell_*.json"))):
    d = json.loads(pathlib.Path(p).read_text(encoding="utf-8"))
    c = d["cell"]
    cells[(c["D"], c["S"], float(c["V"]), c["H"], c["N"])] = d
assert len(cells) == 288


def m(D, S, H="H0", V=0.0, N="bal"):
    return cells[(D, S, V, H, N)]


reps, rois = m("D1", "S2e")["reps"], 27
raw = dict(
    m3_cov_b=m("D1", "S2b")["metrics"]["M3_te0"]["coverage"],
    m1_covered_e=round(m("D1", "S2e")["metrics"]["M1"]["coverage"] * reps * rois, 6),
    m1_covered_b=round(m("D1", "S2b")["metrics"]["M1"]["coverage"] * reps * rois, 6),
    m1_intervals_total=2 * reps * rois,
    d2_h2_s0_slope_reject=m("D2", "S0", "H2")["metrics"]["M3_slope"]["reject"],
    m1_reject_relax_min=min(m("D1", "S2e")["metrics"]["M1"]["reject"], m("D1", "S2b")["metrics"]["M1"]["reject"]),
    m5b_pow_b=m("D1", "S1b")["metrics"]["M5b"]["reject"],
)
raw["m1_covered_total"] = raw["m1_covered_e"] + raw["m1_covered_b"]
rec = dict(source=OUT.name, raw=raw,
           cells_sha256=hashlib.sha256("".join(sorted(hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
                                                      for p in glob.glob(str(OUT / "cells" / "cell_*.json")))).encode()).hexdigest())
(HERE / "s83_values.json").write_text(json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
print(json.dumps(raw, ensure_ascii=False))
