"""Derived ledger value: difference of apparent NDI within one site when only T2 lengthens (Eq. (3), TE weighting only, no new data), with the inputs of Section 2.4 and Brandhofe et al. 2021."""
import json
import math
import pathlib

NU, T2IN, T2EN, K = 0.46, 82.5, 52.0, 113.03 / 102.9


def app(te, ti, tn):
    wi, wn = math.exp(-te / ti), math.exp(-te / tn)
    return NU * wi / (NU * wi + (1 - NU) * wn)


res = dict(inputs=dict(nu=NU, t2_in=T2IN, t2_en=T2EN, ratio=K), by_te={})
for te in (70.0, 100.0):
    c = app(te, T2IN, T2EN)
    res["by_te"][str(int(te))] = dict(control=c, both=app(te, T2IN * K, T2EN * K) - c, extra_only=app(te, T2IN, T2EN * K) - c)
b = [v["both"] for v in res["by_te"].values()]
e = [v["extra_only"] for v in res["by_te"].values()]
res.update(both_abs_min=min(abs(x) for x in b), both_abs_max=max(abs(x) for x in b), extra_abs_min=min(abs(x) for x in e), extra_abs_max=max(abs(x) for x in e))
(pathlib.Path(__file__).resolve().parent / "s74_clinical.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
print(json.dumps(res))
