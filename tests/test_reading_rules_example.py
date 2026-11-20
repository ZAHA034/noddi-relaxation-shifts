"""Tests examples/reading_rules_example.py, pairing each check with a failure control so that every check is able to fail."""
import ast
import importlib.util
import pathlib
import re

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
PATH = ROOT / "examples" / "reading_rules_example.py"
MS = ROOT / "decomp" / "manuscript" / "manuscript.md"


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


EX = _load(PATH, "reading_rules_example")
LIB = _load(ROOT / "decomp" / "stage1" / "lib.py", "stage1_work_lib_for_example")
BANNED = re.compile(r"correct|harmoni[sz]|remov|restor", re.I)


# ---- round trip
def _round_trip_error(forward, n, f):
    x, y = forward(n, f)
    nb, fb = EX.from_log_ratio(x, y)
    return max(np.abs(nb - n).max(), np.abs(fb - f).max())


def _grid():
    rng = np.random.default_rng(1)
    return rng.uniform(0.01, 0.99, 400), rng.uniform(0.01, 0.99, 400)


def test_round_trip_x_y_returns_fractions():
    assert _round_trip_error(EX.to_log_ratio, *_grid()) < 1e-12


def test_failure_control_y_as_logit_fiso_breaks_round_trip():
    def wrong(n, f):
        x, _ = EX.to_log_ratio(n, f)
        return x, np.log(f / (1 - f))
    assert _round_trip_error(wrong, *_grid()) > 0.1


# ---- boundary
def _bounds_ok(eps):
    with np.errstate(divide="ignore", invalid="ignore"):
        x, y = EX.to_log_ratio(np.array([0.0, 1.0, 0.0, 1.0]), np.array([0.0, 1.0, 1.0, 0.0]), eps=eps)
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        return False
    return np.isclose(x[0], np.log(eps / (1 - eps))) and np.isclose(x[1], np.log((1 - eps) / eps))


def test_bounds_0_and_1_clipped_by_eps():
    assert _bounds_ok(1e-4)


def test_failure_control_eps_0_breaks_bounds():
    assert not _bounds_ok(0.0)


# ---- same transform as the analysis code
def test_matches_repository_analysis_function():
    n, f = _grid()
    n, f = np.r_[n, 0.0, 1.0], np.r_[f, 0.0, 1.0]
    c, _ = LIB.coords(n, f, np.zeros_like(n))
    x, y = EX.to_log_ratio(n, f)
    assert LIB.CLIP == 1e-4
    assert np.abs(c[:, 0] - x).max() == 0.0 and np.abs(c[:, 1] - y).max() == 0.0


# ---- model formula and the closed form of proposition 4 (synthetic table · known-answer test · not evidence)
TE, X_G, DR2 = (60.0, 90.0), (0.2, 0.2), (0.012, 0.008)


def _fit(te=TE, x_group=X_G, dr2_group=DR2, with_interaction=True):
    t = EX.synthetic_table(te=te, x_group=x_group, dr2_group=dr2_group)
    x, _ = EX.to_log_ratio(t["ndi"], t["fiso"])
    return EX.fit_site_by_group(x, t["group"], t["site"], with_interaction=with_interaction)


def _matches_closed_form(fit, te=TE, x_group=X_G, dr2_group=DR2):
    want = EX.proposition4_terms(te, x_group, dr2_group)
    return np.isclose(fit.get("group:site", 0.0), want["change"], atol=1e-9) and np.isclose(fit["group"], want["group"], atol=1e-9)


def test_model_has_site_x_group_term():
    assert "group:site" in EX.SITE_BY_GROUP_LINE and "group:site" in EX.SITE_BY_GROUP_REPEATED_LINE
    assert "covariate:site" in EX.SLOPE_BY_SITE_LINE
    assert "(1 | participant)" in EX.SITE_BY_GROUP_REPEATED_LINE


def test_fit_matches_closed_form_of_proposition_4():
    assert _matches_closed_form(_fit())


def test_failure_control_no_interaction_differs_from_closed_form():
    assert not _matches_closed_form(_fit(with_interaction=False))


def test_group_term_keeps_relaxation():
    fit = _fit()
    assert abs(fit["group"]) > 0.1, "분율 차이가 없어도 집단 항은 0 이 아니다(TE × ΔR2 차이)"


def test_equal_TE_gives_zero_site_x_group_term():
    fit = _fit(te=(75.0, 75.0))
    assert abs(fit["group:site"]) < 1e-9 and abs(fit["group"]) > 0.1


def test_fraction_only_difference_gives_zero_site_x_group_term():
    fit = _fit(x_group=(0.2, 0.1), dr2_group=(0.012, 0.012))
    assert abs(fit["group:site"]) < 1e-9 and np.isclose(fit["group"], -0.1)


def test_synthetic_table_fixed_by_seed():
    a, b = EX.synthetic_table(noise_sd=0.02), EX.synthetic_table(noise_sd=0.02)
    assert all(np.array_equal(a[k], b[k]) for k in a)


# ---- the seven «What to do» cells of Table 2 (v16z4: row R3 · differences between acquisitions of groups that share relaxation)
def _table2_rows(text):
    rows, inside = [], False
    for ln in text.split("\n"):
        if ln.startswith("Table: **Table 2.**"):
            inside = True
            continue
        if inside and ln.startswith("| ***"):          # group header rows of Table 2 (2026-11-19, v17h; change record r46/RECORD_gate_2026-11-19.md, a development record outside this repository)
            continue
        if inside and ln.startswith("| ") and not ln.startswith("| Comparison"):
            ln = re.sub(r" \[N-\d+\]", "", ln)        # # ledger tags are not printed (2026-11-20 v17k)
            rows.append([c.strip() for c in ln.strip("| ").split(" | ")])
        elif inside and rows and not ln.startswith("|"):
            break
    return rows


def _rule_mismatches(rules, text):
    rows = _table2_rows(text)
    if len(rows) != len(rules):
        return [f"행 수 {len(rows)} ≠ {len(rules)}"]
    return [r[0] for r, row in zip(rules, rows) if (row[0], row[2]) != (r[0], r[1])]






# ---- words · reading
def test_no_forbidden_words_in_example():
    assert BANNED.findall(PATH.read_text(encoding="utf-8")) == []


def test_failure_control_forbidden_word_check_fires():
    for w in ("corrected", "harmonize", "harmonise", "removes", "restored"):
        assert BANNED.search(w)


def test_reads_no_file_and_uses_only_numpy():
    tree = ast.parse(PATH.read_text(encoding="utf-8"))
    mods = {a.name.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
    mods |= {n.module.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)}
    assert mods <= {"numpy", "__future__"}, mods
    calls = {n.func.id for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    attrs = {n.func.attr for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)}
    assert "open" not in calls and not attrs & {"load", "loadtxt", "genfromtxt", "read_csv", "read_text", "read_bytes"}


def test_run_states_synthetic_and_matches_truth(capsys):
    EX.main()
    out = capsys.readouterr().out
    assert "Synthetic table" in out and "not evidence" in out and "not run here" in out
    assert re.search(r"group term -0\.240 \(closed form -0\.240\), group:site -0\.120 \(closed form -0\.120\)", out)


@pytest.mark.parametrize("word", ["Table 2", "Proposition 4", "known-answer test", "not evidence", "not run"])
def test_docstring_caveat(word):
    assert word in EX.__doc__
