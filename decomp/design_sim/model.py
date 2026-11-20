"""design_sim.model — data generation of the multi-site design simulation (ANALYSIS_PLAN.md Sections 2 · 3), at the level of the signal model without a fitter: y_isr = logit ν_r + θ·g_i + b_ir + TE_s·(R2,en,i − R2,in,i) + γ_sr + δ_sr·ε_isr + η_s·g_i."""
from __future__ import annotations

import itertools
import math

import numpy as np

ACQ = ("OXF1PRI", "OXF2PRI", "OXF3TRI", "NOT2ING-89", "NOT2ING-98", "NOT4GEP", "OXF4GEP", "NOT1ACH")
VENDOR = {"OXF1PRI": "Siemens", "OXF2PRI": "Siemens", "OXF3TRI": "Siemens", "NOT2ING-89": "Philips", "NOT2ING-98": "Philips",
          "NOT4GEP": "GE", "OXF4GEP": "GE", "NOT1ACH": "Philips"}
TE_D1 = {"OXF1PRI": 92.0, "OXF2PRI": 92.0, "OXF3TRI": 96.4, "NOT2ING-89": 89.0, "NOT2ING-98": 98.0, "NOT4GEP": 68.0,
         "OXF4GEP": 69.65, "NOT1ACH": 70.6}                       # card table of the fixed analysis document · ranges at their middle value (same as Supplement S-M9.7)
TE_D2 = {"OXF1PRI": 68.0, "OXF2PRI": 92.0, "OXF3TRI": 98.0, "NOT2ING-89": 69.65, "NOT2ING-98": 89.0, "NOT1ACH": 96.4,
         "NOT4GEP": 70.6, "OXF4GEP": 92.0}                        # short and long TEs for each vendor
TE_C = sum(TE_D1.values()) / len(TE_D1)                           # 84.45625 ms
TR_S = {"OXF1PRI": "3.60", "OXF2PRI": "3.60", "OXF3TRI": "3.60", "NOT2ING-89": "4.00–4.38", "NOT2ING-98": "4.00–4.38",
        "NOT4GEP": "4.90", "OXF4GEP": "8.00", "NOT1ACH": "9.81–10.62"}   # recorded only (TR term 0)
VENDOR_MEAN = {"Siemens": 0.0, "Philips": 0.10, "GE": -0.10}
N_PAT_UNBAL = {"OXF1PRI": 12, "OXF2PRI": 28, "OXF3TRI": 16, "NOT2ING-89": 24, "NOT2ING-98": 20, "NOT4GEP": 20, "OXF4GEP": 28,
               "NOT1ACH": 12}
N_SITE = 40
N_ROI = 27
NU_R = np.linspace(0.45, 0.80, N_ROI)
T2_IN_C, T2_EN_C = 82.5, 52.0
LENGTHEN = 1.0984
SD_B, RHO_ROI = 0.08, 0.5
SD_EPS = 0.031
SD_GAMMA = 0.05
SD_LOGDELTA = 0.1
SD_ETA = 0.02
ETA_SIEMENS = 0.02
CORR_V = 0.8
CV_M5 = 0.10
REPS = 1000
SEED = 20261117

DESIGNS = ("D0", "D1", "D2")
SCENARIOS = ("S0", "S1b", "S1e", "S2b", "S2e", "S2i", "S3b", "S3e")
VS = (0.0, 0.10)
HS = ("H0", "H1", "H2")
NS = ("bal", "unbal")


def delta(assign: str | None, factor: float = LENGTHEN) -> float:
    """Δ = ΔR2,en − ΔR2,in (patients − controls, ms⁻¹) — assignment b both tissue compartments · e extra-neurite · i intra-neurite."""
    if assign is None:
        return 0.0
    fin = factor if assign in ("b", "i") else 1.0
    fen = factor if assign in ("b", "e") else 1.0
    return (1.0 / (T2_EN_C * fen) - 1.0 / T2_EN_C) - (1.0 / (T2_IN_C * fin) - 1.0 / T2_IN_C)


SCEN = {  # name → (assignment the fraction difference is matched to, relaxation assignment)
    "S0": (None, None), "S1b": ("b", None), "S1e": ("e", None), "S2b": (None, "b"), "S2e": (None, "e"), "S2i": (None, "i"),
    "S3b": ("b", "b"), "S3e": ("e", "e"),
}


def theta_of(scen: str) -> float:
    """θ = TE_c · Δ (matched assignment) — the same apparent difference as the relaxation scenario at the central TE (plan Section 3)."""
    m = SCEN[scen][0]
    return 0.0 if m is None else TE_C * delta(m)


def cells() -> list[dict]:
    """The 288 cells — the lexicographic order of (D, S, V, H, N) is the cell number."""
    return [dict(idx=i, D=d, S=s, V=v, H=h, N=n) for i, (d, s, v, h, n) in enumerate(itertools.product(DESIGNS, SCENARIOS, VS, HS, NS))]


def site_te(design: str) -> np.ndarray:
    if design == "D0":
        return np.full(len(ACQ), TE_C)
    tab = TE_D1 if design == "D1" else TE_D2
    return np.array([tab[a] for a in ACQ], float)


def layout(n_design: str):
    """Site number and group (1 patient) of each person — controls first within each site."""
    site, g = [], []
    for s, a in enumerate(ACQ):
        npat = N_SITE // 2 if n_design == "bal" else N_PAT_UNBAL[a]
        site += [s] * N_SITE
        g += [0] * (N_SITE - npat) + [1] * npat
    return np.array(site), np.array(g, float)


def gauss_hermite_raw_theta(theta: float, n: int = 40) -> np.ndarray:
    """For each ROI the true difference on the raw scale E_b[expit(logit ν_r + θ + b) − expit(logit ν_r + b)], b ~ N(0, SD_B²)."""
    x, w = np.polynomial.hermite.hermgauss(n)
    b = math.sqrt(2.0) * SD_B * x
    lv = np.log(NU_R / (1 - NU_R))[:, None]
    ex = lambda z: 1.0 / (1.0 + np.exp(-z))  # noqa: E731
    return (ex(lv + theta + b[None, :]) - ex(lv + b[None, :])) @ w / math.sqrt(math.pi)


def generate(cell: dict, rng: np.random.Generator, *, te_term: bool = True, noise: bool = True, factor: float = LENGTHEN,
             gamma_dev: bool = True, vendor_mean: bool = True, eps: bool = True, delta_noise: bool = True) -> dict:
    """Data of one replicate. Returns: y (people × ROI), site, g, te (site TE of each person), measured r2 (M5a exact · M5b with error), θ, raw-scale θ."""
    site, g = layout(cell["N"])
    n = len(site)
    te_s = site_te(cell["D"])
    te = te_s[site]
    th = theta_of(cell["S"])
    asg = SCEN[cell["S"]][1]
    fin = factor if asg in ("b", "i") else 1.0
    fen = factor if asg in ("b", "e") else 1.0
    t2in = np.where(g == 1, T2_IN_C * fin, T2_IN_C)
    t2en = np.where(g == 1, T2_EN_C * fen, T2_EN_C)
    v = cell["V"]
    if noise and v > 0:
        s = math.sqrt(math.log(1.0 + v * v))
        z1 = rng.standard_normal(n)
        z2 = CORR_V * z1 + math.sqrt(1 - CORR_V ** 2) * rng.standard_normal(n)
        t2in = t2in * np.exp(s * z1 - s * s / 2)
        t2en = t2en * np.exp(s * z2 - s * s / 2)
    r2in, r2en = 1.0 / t2in, 1.0 / t2en
    y = np.tile(np.log(NU_R / (1 - NU_R)), (n, 1)) + th * g[:, None]
    if te_term:
        y = y + (te * (r2en - r2in))[:, None]
    if noise:
        b = SD_B * (math.sqrt(RHO_ROI) * rng.standard_normal((n, 1)) + math.sqrt(1 - RHO_ROI) * rng.standard_normal((n, N_ROI)))
        y = y + b
        gam = np.zeros((len(ACQ), N_ROI))
        if vendor_mean:
            gam += np.array([VENDOR_MEAN[VENDOR[a]] for a in ACQ])[:, None]
        if gamma_dev:
            gam += SD_GAMMA * rng.standard_normal((len(ACQ), N_ROI))
        y = y + gam[site]
        if eps:
            dl = np.exp(SD_LOGDELTA * rng.standard_normal((len(ACQ), N_ROI))) if delta_noise else np.ones((len(ACQ), N_ROI))
            y = y + dl[site] * SD_EPS * rng.standard_normal((n, N_ROI))
        h = cell["H"]
        if h == "H1":
            eta = SD_ETA * rng.standard_normal(len(ACQ))
        elif h == "H2":
            eta = np.array([ETA_SIEMENS if VENDOR[a] == "Siemens" else 0.0 for a in ACQ])
        else:
            eta = np.zeros(len(ACQ))
        y = y + (eta[site] * g)[:, None]
    sm = math.sqrt(math.log(1.0 + CV_M5 ** 2))
    m_in = r2in * np.exp(sm * rng.standard_normal(n) - sm * sm / 2) if noise else r2in
    m_en = r2en * np.exp(sm * rng.standard_normal(n) - sm * sm / 2) if noise else r2en
    return dict(y=y, site=site, g=g, te=te, te_s=te_s, r2=(r2in, r2en), r2_meas=(m_in, m_en), theta=th,
                theta_raw=gauss_hermite_raw_theta(th))
