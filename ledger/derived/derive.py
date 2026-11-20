"""Derived ledger values (2026-11-04): values used in the manuscript that do not come directly from one field of one source file; reads only committed versions of the repository and the fixed copies of the ledger, no images."""
import hashlib, json, math, pathlib, re, subprocess

HERE = pathlib.Path(__file__).resolve().parent
LED = HERE.parent
REPO = LED.parents[4] / "dmri_cov"
DR = LED.parents[1]


def git_show(repo, rel):
    return subprocess.run(["git", "--no-optional-locks", "-C", str(repo), "show", f"HEAD:{rel}"], capture_output=True, check=True).stdout


def src_copy(rel, sha12):
    return (LED / "sources" / rel.replace("../", "__up__/") / sha12).read_bytes()


def theta(t2in, t2en, t2iso):
    return math.degrees(math.atan2(1 / t2in - 1 / t2iso, 1 / t2en - 1 / t2in))


def main():
    inputs = {}
    b = git_show(REPO, "prereg/en/LOCK_B_inputs_EN.md")
    inputs["prereg/en/LOCK_B_inputs_EN.md@HEAD"] = hashlib.sha256(b).hexdigest()
    ln = b.decode("utf-8").splitlines()[76]
    m = re.search(r"T2_in (\d+)±\d+ / (\d+)±\d+ ms, T2_en (\d+)±\d+ / (\d+)±\d+ ms", ln)
    i1, i2, e1, e2 = map(float, m.groups())
    r = git_show(DR, "RECEIPT_2026-11-01b.md")
    inputs["drafts/RECEIPT_2026-11-01b.md@HEAD"] = hashlib.sha256(r).hexdigest()
    m2 = re.search(r"NDI (\d\.\d+) → (\d\.\d+) and FISO (\d\.\d+) → (\d\.\d+)", r.decode("utf-8").splitlines()[17])
    n0, n1, f0, f1 = map(float, m2.groups())

    def card(ndi, fiso):
        return math.log(ndi / (1 - ndi)), math.log(fiso / ((1 - fiso) * ndi))
    (x0, y0), (x1, y1) = card(n0, f0), card(n1, f1)
    gong = {
        "inputs_ms": {"T2_in": [i1, i2], "T2_en": [e1, e2], "T2_iso_assumed": 1000.0},
        "p1_in_mean_en": {"T2_in": i1, "T2_en": (e1 + e2) / 2, "T2_iso": 1000.0, "dir_deg": theta(i1, (e1 + e2) / 2, 1000.0)},
        "mean_both": {"T2_in": (i1 + i2) / 2, "T2_en": (e1 + e2) / 2, "T2_iso": 1000.0, "dir_deg": theta((i1 + i2) / 2, (e1 + e2) / 2, 1000.0)},
        "p1_T2iso2000": {"T2_in": i1, "T2_en": (e1 + e2) / 2, "T2_iso": 2000.0, "dir_deg": theta(i1, (e1 + e2) / 2, 2000.0)},
        "index_shift": {"NDI": [n0, n1], "FISO": [f0, f1], "dx": x1 - x0, "dy": y1 - y0, "dir_deg": math.degrees(math.atan2(y1 - y0, x1 - x0))},
    }
    c = "data/remeasure/explore/cause/cause1b_allunits.json"
    d = "data/remeasure/explore/d5/d5_2_transfer.json"
    cb, db = src_copy(c, "6f69f3710e88"), src_copy(d, "4375d393a348")
    inputs[c] = hashlib.sha256(cb).hexdigest()
    inputs[d] = hashlib.sha256(db).hexdigest()
    cj, dj = json.loads(cb), json.loads(db)
    sd_ = dj["sigma_d_median"]["x"]
    retest = {"x_mean": cj["sets"]["all"]["x"]["mean"], "x_sd": cj["sets"]["all"]["x"]["sd"], "sigma_d_x": sd_,
              "mean_over_sigma_d": cj["sets"]["all"]["x"]["mean"] / sd_, "sd_over_sigma_d": cj["sets"]["all"]["x"]["sd"] / sd_}
    rb = "data/remeasure/explore/rv/rv_b_family.json"
    rbb = src_copy(rb, "747df945ec53")
    inputs[rb] = hashlib.sha256(rbb).hexdigest()
    obs = json.loads(rbb)["obs"]
    for key in ("p1_in_mean_en", "mean_both"):                 # Observed TE slope (s⁻¹) ÷ logarithmic slope of the same literature T2
        g = gong[key]
        g["slope_ratio_x"] = obs["bx"] / ((1 / g["T2_en"] - 1 / g["T2_in"]) * 1000)
        g["slope_ratio_y"] = obs["by"] / ((1 / g["T2_in"] - 1 / g["T2_iso"]) * 1000)
        g["family_T2_en"] = 1 / (obs["bx"] / 1000 + 1 / g["T2_in"])    # Set of solutions that fit the two observed slopes (signal model) — T2_en · T2_iso at this T2_in
        g["family_T2_iso"] = 1 / (1 / g["T2_in"] - obs["by"] / 1000)
    out = {"what": "장부 도출 값(derive.py)", "inputs_sha256": inputs, "gong": gong, "retest": retest}
    (HERE / "derived.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: v.get("dir_deg") for k, v in gong.items() if isinstance(v, dict) and "dir_deg" in v}), retest["mean_over_sigma_d"], retest["sd_over_sigma_d"])


if __name__ == "__main__":
    main()
