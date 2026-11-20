import importlib.util, numpy as np
spec = importlib.util.spec_from_file_location("d5", "<WORK>/stage1_work/drafts/x14/D5/d5_analysis.py")
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)
rng = np.random.default_rng(1)
# Structure similar to the real one: 7 devices · sizes 30/20/20/20/10/10/9 · 20 subjects
sizes = [30, 20, 20, 20, 10, 10, 9]; g = np.concatenate([[f"d{i}"] * s for i, s in enumerate(sizes)]); n = len(g)
subj = np.concatenate([np.arange(s) % 20 for s in sizes]).astype(str)
p = 4
cut = A.null_cutoff(g, subj, p, 0.99, 100, 5); cut_sub = A.null_cutoff(g, subj, p, 0.99, 100, 6, sub=True)
from scipy.stats import chi2
print("cut", round(cut, 2), "cut_sub", round(cut_sub, 2), "chi2_99", round(chi2.ppf(0.99, p), 2))
# Null-hypothesis data: normal · stable outlier count and simple outlier count
fs, ss = [], []
for s in range(20):
    r2 = np.random.default_rng(200 + s); V = r2.normal(0, 1, (n, p)) * r2.uniform(0.5, 3, len(sizes))[np.searchsorted(sorted(set(g)), g)][:, None]
    d2, fl = A.mcd_flags(V, g, cut); fq, _ = A.stability(V, g, subj, 50, 300 + s, cut_sub); fs.append(int(fl.sum())); ss.append(int(((fq >= 0.8) & fl).sum()))
print("null flags mean", np.mean(fs), "(expected", 0.01 * n, ") stable mean", np.mean(ss), "max", max(ss))
# Planted outlier set (+6 SD on all axes)
V = rng.normal(0, 1, (n, p)); V[[3, 60, 110]] += 6
d2, fl = A.mcd_flags(V, g, cut); fq, _ = A.stability(V, g, subj, 50, 7, cut_sub)
st = (fq >= 0.8) & fl
print("planted found", st[[3, 60, 110]].tolist(), "other stable", int(np.delete(st, [3, 60, 110]).sum()))
assert st[[3, 60, 110]].all() and np.mean(ss) <= 1.5
print("ALL OK")
