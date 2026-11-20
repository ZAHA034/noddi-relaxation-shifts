"""S73 §5.3 and §5.4 readings from the stored results (no new computation)."""
from __future__ import annotations

import json
import os
import pathlib
import sys

REPO = pathlib.Path(os.environ.get("S73_REPO") or pathlib.Path(__file__).resolve().parents[3])
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import rconfig as R  # noqa: E402,F401


def category(r: float) -> str:
    return "설명하지 못함" if r <= 0.25 else "일부 설명" if r < 0.75 else "대부분 설명할 수 있음"


def row_5_4(beta: dict, ci: dict) -> str | None:
    """Conditions of the implementation note — first two rows: β > 0 for all three participants and no interval contains 0 · third row: the interval of at least one participant contains 0, and no participant has β < 0 with an interval that excludes 0 · otherwise: None."""
    pos = all(beta[s] > 0 and ci[s][0] > 0 for s in beta)
    neg_sig = any(beta[s] < 0 and ci[s][1] < 0 for s in beta)
    inc0 = any(ci[s][0] <= 0 <= ci[s][1] for s in beta)
    if pos:
        return "first_or_second"
    if inc0 and not neg_sig:
        return "third"
    return None


def main(d: str) -> dict:
    out = pathlib.Path(d)
    noise = json.loads((out / "s73_noise.json").read_text(encoding="utf-8"))
    align = json.loads((out / "s73_align.json").read_text(encoding="utf-8"))
    res = dict(definition="decomp/plans/S73_POST_HOC_PLAN_2026-11-16.md §2.5 · §5.3 · §5.4 · S73_AMENDMENT_2026-11-16.md",
               section2_stopped=align.get("stopped"))
    paths = align.get("paths", {})
    full = {p: isinstance(v.get("env"), dict) and "tau_pooled_x" in v["env"] for p, v in paths.items()}
    res["paths_with_e1"] = [p for p, ok in full.items() if ok]
    if not paths or not full.get("A"):
        res["note"] = "경로가 계산되지 않음 — 판독하지 않는다"
        (out / "s73_reading.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
        return res
    # §5.3
    tau_obs_by = {p: paths[p]["env"]["tau_pooled_x"] for p in res["paths_with_e1"]}
    p_min = min(tau_obs_by, key=tau_obs_by.get)
    tau_sim = {L: v["pooled"]["tau"] for L, v in noise["by_L"].items()}
    Lt = max(tau_sim, key=tau_sim.get)
    res["tau"] = dict(obs_by_path=tau_obs_by, obs=tau_obs_by[p_min], obs_path=p_min, sim_by_L=tau_sim, L_max=Lt,
                      ratio=tau_sim[Lt] / tau_obs_by[p_min], category=category(tau_sim[Lt] / tau_obs_by[p_min]))
    # β of §5.4 (path B, otherwise C)
    pb = "B" if full.get("B") else "C" if full.get("C") else None
    res["beta_path"] = pb
    subs = sorted(noise["by_L"]["1"]["by_participant"])
    if pb:
        bp = paths[pb]["env"]["by_participant"]
        beta = {s: bp[s]["beta"] for s in subs}
        ci = {s: bp[s]["beta_ci95"] for s in subs}
        res["beta"] = dict(path=pb, beta=beta, ci95=ci, row=row_5_4(beta, ci))
        r = res["beta"]["row"]
        res["row_5_4"] = (None if r is None else "third" if r == "third" else
                          ("first" if res["tau"]["category"] == "설명하지 못함" else "second"))
        res["row_5_4_part_or_most"] = (None if res["row_5_4"] != "second" else "most" if res["tau"]["category"] == "대부분 설명할 수 있음" else "part")
        # β, simulated against observed (same bounds)
        def bsim(L, s):
            return noise["by_L"][L]["by_participant"][s]["beta"]
        res["beta_sim_vs_obs"] = {s: dict(largest=max(bsim(L, s) for L in noise["by_L"] if bsim(L, s) == bsim(L, s)) / beta[s],
                                          at_L_of_largest_tau=bsim(Lt, s) / beta[s]) for s in subs}
        res["odi_sim_vs_obs"] = {s: max(noise["by_L"][L]["by_participant"][s]["odi_range_ratio"] for L in noise["by_L"])
                                 / paths[pb]["odi_range_ratio"][s] for s in subs}
    # §2.5
    po = "B" if full.get("B") else "C" if full.get("C") else None
    s25 = dict(path=po)
    if po:
        ra, rp = paths["A"]["odi_range_ratio"], paths[po]["odi_range_ratio"]
        q = {s: rp[s] / ra[s] for s in subs}
        s25["odi_ratio_path_over_A"] = q
        s25["odi_path_range"] = [min(rp.values()), max(rp.values())]
        s25["odi_ending"] = ("registration accounted for most of it" if all(v <= 0.5 for v in q.values()) else
                             "separate registration did not account for it" if all(v >= 0.8 for v in q.values()) else
                             "accounted for part of it")
        dA, dP = paths["A"]["env"]["dte_star_x"], paths[po]["env"]["dte_star_x"]
        s25["dte_star"] = dict(A=dA, path=dP)
        e4a, e4p = paths["A"]["env"]["e4_wong_62_132_x_minus_raw"], paths[po]["env"]["e4_wong_62_132_x_minus_raw"]
        s25["e4"] = dict(A=e4a, path=e4p, sign_changed=(e4a["diff"] > 0) != (e4p["diff"] > 0))
        s25["psi"] = dict(A=paths["A"]["psi"]["overall"]["psi_median"], path=paths[po]["psi"]["overall"]["psi_median"],
                          A_ci90=paths["A"]["psi"]["ci90"], path_ci90=paths[po]["psi"]["ci90"])
    lin = {p: paths[p]["linear"] for p in paths}
    s25["linear"] = lin
    res["section_2_5"] = s25
    (out / "s73_reading.json").write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    return res


if __name__ == "__main__":
    print(json.dumps(main(sys.argv[1]), ensure_ascii=False, indent=1))
