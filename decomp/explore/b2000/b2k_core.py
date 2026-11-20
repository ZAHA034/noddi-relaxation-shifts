"""Shared mask and metric code for the b2000 noise floor, used by the measurement (b2k_read) and the known-answer test (b2k_sim)."""
from __future__ import annotations

import math

import numpy as np
from scipy import ndimage as ndi
from scipy import special

B0_MAX = 50.0                 # b ≤ 50 = b0 (same as rconfig.B0_MAX_PREP)
B2K_CENTER, B2K_HALF = 2000.0, 100.0
BOX = dict(R=(-0.25, 0.25), A=(-0.35, 0.35), S=(-0.20, 0.40))
SMOOTH_MM = 15.0
Q_THR = 1.5
MIN_VOX = 50


def shells(bvals) -> tuple[np.ndarray, np.ndarray]:
    bv = np.asarray(bvals, float).ravel()
    return np.flatnonzero(bv <= B0_MAX), np.flatnonzero(np.abs(bv - B2K_CENTER) <= B2K_HALF)


def ras_grid(shape: tuple[int, int, int], affine: np.ndarray) -> np.ndarray:
    """RAS coordinates (mm) of the voxel centers — (nx, ny, nz, 3)."""
    ijk = np.stack(np.meshgrid(*[np.arange(n, dtype=np.float64) for n in shape], indexing="ij"), -1)
    return ijk @ affine[:3, :3].T + affine[:3, 3]


def otsu(x: np.ndarray, nbins: int = 256) -> float:
    """Otsu threshold (histogram between the 0.5 and 99.5 percentiles)."""
    x = np.asarray(x, np.float64).ravel()
    lo, hi = np.percentile(x, [0.5, 99.5])
    h, e = np.histogram(np.clip(x, lo, hi), bins=nbins, range=(lo, hi))
    c = (e[:-1] + e[1:]) / 2
    w0 = np.cumsum(h).astype(np.float64)
    w1 = w0[-1] - w0
    s0 = np.cumsum(h * c)
    m0 = s0 / np.maximum(w0, 1)
    m1 = (s0[-1] - s0) / np.maximum(w1, 1)
    return float(e[1 + int(np.argmax(w0 * w1 * (m0 - m1) ** 2))])


def head_mask(B0: np.ndarray) -> np.ndarray:
    """[fix before reading §27-1] threshold = Otsu(B0). The first definition (0.2 × the 99th percentile) took the background floor for the head when L was large and SNR low (known-answer test)."""
    m = B0 > otsu(B0)
    lab, n = ndi.label(m, structure=np.ones((3, 3, 3), bool))
    if n == 0:
        return m
    sizes = ndi.sum(m, lab, np.arange(1, n + 1))
    return ndi.binary_fill_holes(lab == 1 + int(np.argmax(sizes)))


def central_box(H: np.ndarray, X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    P = X[H]
    c = P.mean(0)
    L = np.percentile(P, 98, axis=0) - np.percentile(P, 2, axis=0)
    d = X - c
    C = np.ones(H.shape, bool)
    for ax, key in enumerate("RAS"):
        lo, hi = BOX[key]
        C &= (d[..., ax] >= lo * L[ax]) & (d[..., ax] <= hi * L[ax])
    return C, c, L


def local_contrast(B0: np.ndarray, H: np.ndarray, zooms) -> np.ndarray:
    """Q = B0 / Ŝ · Ŝ = G(B0·H)/G(H) — divides out the coil-sensitivity bias with a normalized Gaussian (σ 15 mm)."""
    sig = [SMOOTH_MM / float(z) for z in zooms[:3]]
    num = ndi.gaussian_filter(np.where(H, B0, 0.0), sig, mode="constant")
    den = ndi.gaussian_filter(H.astype(np.float64), sig, mode="constant")
    S = np.divide(num, den, out=np.zeros_like(num), where=den > 1e-6)
    return np.divide(B0, S, out=np.zeros_like(S), where=S > 0)


def ventricle_mask(B0: np.ndarray, affine: np.ndarray, zooms) -> tuple[np.ndarray, dict]:
    B0 = np.asarray(B0, np.float64)
    H = head_mask(B0)
    X = ras_grid(B0.shape, np.asarray(affine, float))
    C, c, L = central_box(H, X)
    Q = local_contrast(B0, H, zooms)
    K = H & C & (Q >= Q_THR)
    V = ndi.binary_erosion(K, structure=ndi.generate_binary_structure(3, 1), iterations=1)
    info = dict(n_H=int(H.sum()), n_C=int((H & C).sum()), n_K=int(K.sum()), n_V=int(V.sum()),
                centroid_mm=[round(float(v), 2) for v in c], extent_mm=[round(float(v), 2) for v in L])
    return V, info


def floor_metrics(b2k: np.ndarray, V: np.ndarray, b0: np.ndarray | None = None) -> dict:
    """Main metric R · secondary metric S1. b2k = (…, N2) b2000 magnitude values · b0 = (…, N0) — voxels that are 0 or below in any volume and voxels with SD 0 are dropped."""
    X = np.asarray(b2k, np.float64)[V]
    bad = (X <= 0).any(1)
    if b0 is not None:
        bad |= (np.asarray(b0, np.float64)[V] <= 0).any(1)
    sd = X.std(1, ddof=1)
    good = ~bad & (sd > 0)
    out = dict(n_V=int(V.sum()), n_used=int(good.sum()), n_nonpos=int(bad.sum()), n_sd0=int((~bad & (sd <= 0)).sum()),
               N2=int(X.shape[1]) if X.ndim == 2 else 0)
    if out["n_V"] < MIN_VOX or good.sum() < MIN_VOX:
        out.update(R=None, status="마스크 부족")
        return out
    Xg = X[good]
    r = Xg.mean(1) / sd[good]
    q = np.percentile(r, [25, 50, 75])
    s1 = Xg.mean(0) / Xg.std(0, ddof=1)
    out.update(R=float(q[1]), r_q25=float(q[0]), r_q75=float(q[2]), S1=float(np.median(s1)),
               mean_M=float(Xg.mean()), mean_M2=float((Xg ** 2).mean()), status="값")
    return out


# ── Theory · adjustment ──────────────────────────────────────────────────────────────
def chi_ratio(L: float) -> float:
    """Mean/SD of the central chi (2L degrees of freedom · signal 0) — population values (for checking)."""
    k = 2.0 * L
    m = math.sqrt(2.0) * math.exp(special.gammaln((k + 1) / 2) - special.gammaln(k / 2))
    return m / math.sqrt(k - m * m)


def ncchi(s: np.ndarray, sigma: np.ndarray | float, L: float, rng: np.random.Generator) -> np.ndarray:
    """Noncentral chi magnitude — M = √((s + σz)² + σ² · Gamma((2L − 1)/2, 2)) · for L = 0.5 it is |s + σz| (real output)."""
    s = np.asarray(s, np.float64)
    z = rng.standard_normal(s.shape)
    re = s + sigma * z
    if L <= 0.5 + 1e-12:
        return np.abs(re)
    g = rng.gamma((2.0 * L - 1.0) / 2.0, 2.0, size=s.shape)
    return np.sqrt(re * re + (np.asarray(sigma) ** 2) * g)


def l_eff(R: float, Ls, Rmed) -> float | None:
    """Inverts the calibration curve (L → median R · monotone) by linear interpolation in log L — outside the range a boundary mark (±inf) is returned."""
    Ls, Rmed = np.asarray(Ls, float), np.asarray(Rmed, float)
    if R is None or not np.all(np.diff(Rmed) > 0):
        return None
    if R < Rmed[0]:
        return -math.inf
    if R > Rmed[-1]:
        return math.inf
    return float(np.exp(np.interp(R, Rmed, np.log(Ls))))
