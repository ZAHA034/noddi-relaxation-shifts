import sys, importlib.util, numpy as np
spec = importlib.util.spec_from_file_location("d1", "<WORK>/stage1_work/drafts/x14/D1/d1_c3.py")
D = importlib.util.module_from_spec(spec); spec.loader.exec_module(D)
J = D.J
rng = np.random.default_rng(0)
# 1) Rec — fake draw_fn: a different mean for each combination j · k
n, m, nk = 5, 2, 3
calls = []
def fake(X, cond, n_keep, gs):
    j, kk = divmod(len(calls), nk); calls.append((j, kk))
    T = rng.normal(0, 1, (n, n_keep, 6)).astype(np.float32) + (j * 10 + kk)
    ok = np.ones(n, bool); ok[0] = (j == 0)
    T[~ok] = np.nan
    return T, ok, np.full(n, 0.9)
rec = D.Rec(fake, n, m, nk)
for j in range(m):
    for kk in range(nk):
        out = rec(None, None, 400, None)
assert rec.c == m * nk
assert np.allclose(np.nanmean(rec.mean[1:, 1, 2, :]), 12, atol=0.2), rec.mean[1:, 1, 2]
assert np.isnan(rec.pit[0, 1]).all() and not rec.ok[0, 1].any() and rec.ok[0, 0].all()
assert (rec.pit[1:, 1] < 0.01).all() and np.allclose(rec.pit[1:, 0, 0], 0.5, atol=0.1)
print("1 Rec ok")
# 2) pick_rule — if the PIT of combination 1 is farther from 0.5, the z of that combination
pit = np.full((3, 2, 9, 3), 0.5); pit[:, 1] = 0.9; pit[0, 0] = 0.01
mean = np.zeros((3, 2, 9, 3)); mean[:, 0] = -1; mean[:, 1] = 2; sd = np.ones((3, 2, 9, 3))
P, Z, j = D.pick_rule(pit, mean, sd, np.zeros((3, 3)))
assert (P[0] == 0.01).all() and (Z[0] == -1).all() and (P[1:] == 0.9).all() and (Z[1:] == 2).all()
pit[2, 0] = 0.1; pit[2, 1] = 0.9   # Same distance → earlier combination
P, Z, j = D.pick_rule(pit, mean, sd, np.zeros((3, 3)))
assert (P[2] == 0.1).all() and (Z[2] == -1).all()
print("2 pick_rule ok")
# 3) decomp — 4 ROIs × 50 subjects · SD within ROI 1 · SD of ROI means 3
rois = np.repeat(np.arange(4), 50); V = rng.normal(0, 1, (200, 9)) + np.array([0, 3, -3, 6])[rois][:, None]
a, b, e = D.decomp(V, rois.astype(str), np.ones(200, bool))
assert np.allclose(a, 1, atol=0.15), a
assert np.allclose(b, np.std([0, 3, -3, 6], ddof=1), atol=0.4), b
assert np.allclose(e, a / np.sqrt(50))
print("3 decomp ok")
# 4) kstats
Pk = rng.uniform(0, 1, (1000, 9)); Zk = rng.normal(0, 1, (1000, 9))
s = D.kstats(Pk, Zk, np.ones(1000, bool))
assert np.allclose(s["pit_mean"], 0.5, atol=0.05) and np.allclose(s["tail"], 0.10, atol=0.03) and np.allclose(s["z_sd"], 1, atol=0.08)
print("4 kstats ok")
# 5) compare — stored values made from the same data are equal · changing the coverage of one ROI unit breaks the match
nn = 60
u = dict(cov=rng.uniform(size=(nn, 9, 3)) < 0.8, ok=rng.uniform(size=(nn, 9)) < 0.97, flag=rng.uniform(size=nn) < 0.1,
         subjects=np.repeat(np.array(["s1", "s2", "s3", "s4"]), 15), main_excluded=np.zeros(nn, bool),
         combos=["32|32", "32|64"], flag_by_combo=rng.uniform(size=(nn, 2)) < 0.2, ok_by_combo=rng.uniform(size=(nn, 2)) < 0.95)
us = J.select(u, ~u["main_excluded"])
c4 = J.cell_verdict(us["cov"][:, :, 2], us["ok"], us["flag"], us["subjects"], 0.8, np.random.default_rng(1), 0.2)
c4["by_coil_combo"] = J.coil_counts(us)
import json; c4 = json.loads(json.dumps(c4))
r = D.compare(u, c4, 2); assert r["match"] and r["mismatch_units"] == 0 and r["mismatch_coil"] == 0, r
u2 = dict(u); u2["cov"] = u["cov"].copy(); i = int(np.flatnonzero(~u["flag"] & u["ok"].all(1))[0]); u2["cov"][i, 0, 2] = ~u2["cov"][i, 0, 2]
r2 = D.compare(u2, c4, 2); assert not r2["match"] and r2["mismatch_units"] == 1, r2
u3 = dict(u); u3["flag_by_combo"] = u["flag_by_combo"].copy(); u3["flag_by_combo"][0, 0] = ~u3["flag_by_combo"][0, 0]
r3 = D.compare(u3, c4, 2); assert not r3["match"] and r3["mismatch_coil"] == 1, r3
r4 = D.compare(dict(cov=np.zeros((0, 9, 3), bool)), {"verdict": "자료 없음"}, 2); assert r4["match"]
print("5 compare ok")
print("combos_of", D.combos_of("OXF1PRI~OXF2PRI"), D.combos_of("OXF1PRI~NOT1ACH"), len(D.combos_of("NOT4GEP~NOT1ACH")))
print("needed types", D.needed_types())
print("ALL OK")
