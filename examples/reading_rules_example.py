"""Applying the reading rules of Table 2 to a table of ROI values: a short worked example.

This is an example of applying the reading rules of Table 2 to a table of ROI values; it is not a method and does not
separate relaxation from tissue. The site × group term contains the relaxation-driven residual of Proposition 4 together
with any true site-dependent group difference; telling them apart requires compartmental relaxation measured in each
group at each site. The group term still includes relaxation: at each site the group difference contains TE times the
group difference in ΔR2. A site × group term near zero does not show that the groups share relaxation, because the
residual vanishes at equal TE.

These statements hold under the assumptions of Proposition 4: the groups share T1 and proton density, and there is no
multiplicative site term.

The fit to the synthetic table below is a known-answer test of this code, not evidence: the synthetic table is generated
from the data model of Proposition 4. Every value in it comes from a fixed seed; none comes from a study, and none is used
in the article.

The model lines are written for R or for statsmodels. Only the fixed-effects line of the last rule is fitted here (numpy
least squares), so that the example needs nothing beyond numpy; the line with a participant random intercept is a string
and is not run.

Run:  python examples/reading_rules_example.py
"""
from __future__ import annotations

import numpy as np

# Table 2, row by row: comparison, the «What to do» cell as printed, and a model line where the cell can be written as one.
READING_RULES = (
    ("Association with a covariate",
     "Use logit NDI as the outcome; rank people only if their relaxation varies little",
     "x ~ covariate + covariates  (within one protocol)"),
    ("Patient–control difference",
     "Report d/TE; state the assumption or measure relaxation",
     None),
    ("Levels (absolute NDI)",
     "Compare within-acquisition slopes of logit NDI",
     "x ~ covariate + site + covariate:site + covariates"),
    ("Group difference, shared relaxation",
     "Compare logit NDI differences; state the assumption",
     "x ~ group + site + covariates  (groups sharing relaxation)"),
    ("Group difference after additive site adjustment",
     "Measure relaxation (for example, a second TE), or fit group × TE and read the difference extrapolated to zero TE "
     "(shared T1 and proton density; SE about 7 times larger); a site × group term alone shows the between-site residual, "
     "not the fraction difference",
     "x ~ group + site + group:site + covariates"),
    ("FISO",
     "Compare at equal TR and TE; site term across sites; measure relaxation before reading free water",
     None),
    ("ODI",
     "Treat scanner model and protocol as a factor",
     "odi ~ group + scanner_protocol + covariates"),
)
SITE_BY_GROUP_LINE = next(r[2] for r in READING_RULES if r[0] == "Group difference after additive site adjustment")
SLOPE_BY_SITE_LINE = READING_RULES[2][2]
# With repeated measures of a participant. A string: not run here.
SITE_BY_GROUP_REPEATED_LINE = SITE_BY_GROUP_LINE + " + (1 | participant)"


def to_log_ratio(ndi, fiso, eps=1e-4):
    """x = logit NDI and y = logit FISO − ln NDI, with both fractions kept at least eps from the bounds (Section 3.2).

    Apply it after ROI averaging: the fractions in the table are ROI means (Section 3.2). The translation holds for each
    tissue state, so it is approximate for fractions averaged over heterogeneous tissue before the transform.
    """
    n = np.clip(np.asarray(ndi, dtype=float), eps, 1 - eps)
    f = np.clip(np.asarray(fiso, dtype=float), eps, 1 - eps)
    return np.log(n / (1 - n)), np.log(f / (1 - f)) - np.log(n)


def from_log_ratio(x, y):
    """The fractions (NDI, FISO) back from x and y; used by the round-trip test."""
    n = np.exp(x) / (1 + np.exp(x))
    logit_f = np.asarray(y, dtype=float) + np.log(n)
    return n, np.exp(logit_f) / (1 + np.exp(logit_f))


def synthetic_table(te=(60.0, 90.0), x_group=(0.2, 0.2), dr2_group=(0.012, 0.008), tau=(0.0, 0.3), n_per_cell=30,
                    noise_sd=0.0, seed=20261108):
    """A table of ROI values generated from the data model of Proposition 4, x = x_g + TE_s ΔR2_g + τ_s (+ noise).

    Index 0 of x_group and dr2_group is the control group, index 1 the patient group; te and tau give the two sites. ΔR2 is
    in 1/ms. Returns a dict of equal-length arrays: site (0 or 1), group (0 or 1), te, ndi, fiso.
    """
    rng = np.random.default_rng(seed)
    site, group, te_col, x = [], [], [], []
    for s, (te_s, tau_s) in enumerate(zip(te, tau)):
        for g, (x_g, dr2_g) in enumerate(zip(x_group, dr2_group)):
            site += [s] * n_per_cell
            group += [g] * n_per_cell
            te_col += [te_s] * n_per_cell
            x += list(x_g + te_s * dr2_g + tau_s + rng.normal(0.0, noise_sd, n_per_cell))
    x = np.array(x)
    return dict(site=np.array(site), group=np.array(group), te=np.array(te_col), ndi=np.exp(x) / (1 + np.exp(x)),
                fiso=rng.uniform(0.05, 0.20, x.size))


def fit_site_by_group(x, group, site, with_interaction=True):
    """Least-squares coefficients of the fixed-effects line of the site × group rule, with site 0 as the reference site.

    Returns a dict: «group» is the group difference at site 0, «group:site» its change from site 0 to site 1.
    """
    names = ["intercept", "group", "site"] + (["group:site"] if with_interaction else [])
    columns = [np.ones_like(x), group, site] + ([group * site] if with_interaction else [])
    beta = np.linalg.lstsq(np.column_stack(columns).astype(float), x, rcond=None)[0]
    return dict(zip(names, beta))


def proposition4_terms(te, x_group, dr2_group):
    """The closed form of Proposition 4 for the two sites: the group term at site 0 and the between-site change."""
    ddr2 = dr2_group[1] - dr2_group[0]
    return dict(group=(x_group[1] - x_group[0]) + te[0] * ddr2, change=(te[1] - te[0]) * ddr2)


def main():
    print("Reading rules of Table 2 and a model line where the rule can be written as one:")
    for comparison, what_to_do, line in READING_RULES:
        print(f"- {comparison}\n    what to do: {what_to_do}\n    model line: {line or 'none: state or measure, as the cell says'}")
    print(f"With repeated measures (string only, not run here): {SITE_BY_GROUP_REPEATED_LINE}")
    print("\nSynthetic table, generated from Proposition 4 (a known-answer test of this code, not evidence):")
    te, x_group, dr2_group = (60.0, 90.0), (0.2, 0.2), (0.012, 0.008)
    table = synthetic_table(te=te, x_group=x_group, dr2_group=dr2_group)
    x, y = to_log_ratio(table["ndi"], table["fiso"])
    fit = fit_site_by_group(x, table["group"], table["site"])
    want = proposition4_terms(te, x_group, dr2_group)
    print(f"  no fraction difference, relaxation differs: group term {fit['group']:+.3f} (closed form {want['group']:+.3f}), "
          f"group:site {fit['group:site']:+.3f} (closed form {want['change']:+.3f})")
    print(f"  y computed for {y.size} rows; synthetic values, not for any report")


if __name__ == "__main__":
    main()
