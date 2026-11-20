"""Recomputes the two «TR on FISO» values of Section 2.4 of an earlier manuscript version (numbering of that version) with repository code, on CPU, writing nothing to the repository."""
import itertools
import json
import sys

import numpy as np

sys.dont_write_bytecode = True
repo = sys.argv[1]
sys.path[:0] = [f"{repo}/decomp/stage1", f"{repo}/decomp/s4"]
import config as C  # noqa: E402
import s4cfg as S  # noqa: E402

TR1, TR2, F_TRUE = 3600.0, 9800.0, 0.10
grid = dict(c=np.linspace(*S.T1C_RANGE, 26), nu=np.linspace(*S.BOX["NDI"], 26),
            te=np.linspace(60.0, 110.0, 11), t2in=np.linspace(*S.T2IN_RANGE, 6), r=np.linspace(*S.R_RANGE, 9))


def t1s(c, ratio):
    return (c * np.sqrt(ratio), c / np.sqrt(ratio), S.T1ISO)


out = {}
for scen, ratio in S.SCEN.items():
    dl, draw = [], []
    for c, nu in itertools.product(grid["c"], grid["nu"]):
        t1 = t1s(c, ratio)
        for te, t2in, r in itertools.product(grid["te"], grid["t2in"], grid["r"]):
            t2 = (t2in, t2in / r, S.T2ISO)
            a1, a2 = C.weights(te, t2, TR1, t1), C.weights(te, t2, TR2, t1)
            dl.append(np.log(a2[2] / a1[2]) - np.log((nu * a2[0] + (1 - nu) * a2[1]) / (nu * a1[0] + (1 - nu) * a1[1])))
        w1, w2 = 1 - np.exp(-TR1 / np.asarray(t1)), 1 - np.exp(-TR2 / np.asarray(t1))
        fa = [F_TRUE * w[2] / (F_TRUE * w[2] + (1 - F_TRUE) * (nu * w[0] + (1 - nu) * w[1])) for w in (w1, w2)]
        draw.append(fa[1] - fa[0])
    out[scen] = dict(t1_ratio=ratio, dlogit_min=float(min(dl)), dlogit_max=float(max(dl)),
                     draw_min=float(min(draw)), draw_max=float(max(draw)))
iso = float(np.log((1 - np.exp(-TR2 / S.T1ISO)) / (1 - np.exp(-TR1 / S.T1ISO))))
res = dict(TR=[TR1, TR2], f_true=F_TRUE, grid={k: [float(v[0]), float(v[-1]), len(v)] for k, v in grid.items()},
           by_scenario=out, all=dict(dlogit=[min(v["dlogit_min"] for v in out.values()), max(v["dlogit_max"] for v in out.values())],
                                     draw=[min(v["draw_min"] for v in out.values()), max(v["draw_max"] for v in out.values())]),
           free_water_term_upper_bound=iso)
print(json.dumps(res, indent=1))
