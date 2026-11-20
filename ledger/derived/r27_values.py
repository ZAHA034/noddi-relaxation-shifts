"""R27 (v16z3): the two ratios of the retest-SD sentence of Section 4.1 (median retest SD of c10 ÷ median retest SD of OXF1PRI, for x and y); two divisions, no new data."""
import hashlib
import json
import pathlib
import sys
from decimal import Decimal, ROUND_HALF_UP

sys.dont_write_bytecode = True
REPO = pathlib.Path(sys.argv[1]).resolve()
HERE = pathlib.Path(__file__).resolve().parent
SRC = {
    "judgements": "decomp/results/remeasure/judgements.json",
    "s73_fitter": "data/remeasure/explore/s73/s73_fitter.json",
    "s1a": "<project-notes>/x15/env/out/s1a.json",
}


def real(rel):
    return (REPO / rel).resolve()


def load(rel):
    b = real(rel).read_bytes()
    return json.loads(b.decode("utf-8")), hashlib.sha256(b).hexdigest()


def rnd(x, q):
    return format(Decimal(repr(x)).quantize(Decimal(q), rounding=ROUND_HALF_UP), "f")


J, hJ = load(SRC["judgements"])
F, hF = load(SRC["s73_fitter"])
A, hA = load(SRC["s1a"])
c10_x, c10_y = J["Tong"]["main"]["s_xy_unit_median"]
onh_x, onh_y = F["retest_sd"]["x"], F["retest_sd"]["y"]
if A["s_x_median"] != onh_x:
    raise SystemExit(f"대조 (1) 실패: s1a {A['s_x_median']} ≠ s73_fitter {onh_x}")
got = [rnd(v, "0.001") for v in (c10_x, c10_y, onh_x, onh_y)]
if got != ["0.111", "0.064", "0.031", "0.084"]:
    raise SystemExit(f"대조 (2) 실패: {got}")
ratio_x, ratio_y = c10_x / onh_x, c10_y / onh_y
printed = dict(ratio_x=rnd(ratio_x, "0.1"), ratio_y=rnd(ratio_y, "0.01"))
if [rnd(c10_x / onh_y, "0.1"), rnd(c10_y / onh_x, "0.01")] == [printed["ratio_x"], printed["ratio_y"]]:
    raise SystemExit("대조 (3) 실패: 입력을 바꿔도 인쇄 값이 같음")
out = dict(
    what="4.1 재검 SD 문장의 두 비(c10 ÷ OXF1PRI) · R27(v16z3) · 저자 결정(1i)",
    inputs=dict(c10_x=c10_x, c10_y=c10_y, onh_x=onh_x, onh_y=onh_y),
    ratio_x=ratio_x, ratio_y=ratio_y, printed=printed,
    sources={k: dict(path=v, sha256=h) for (k, v), h in zip(SRC.items(), (hJ, hF, hA))},
    controls=dict(s1a_equals_s73_x=True, printed_sds=got, swapped_inputs_change_printed=True),
)
(HERE / "r27_values.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
print(json.dumps(dict(printed=printed, sds=got), ensure_ascii=False))
