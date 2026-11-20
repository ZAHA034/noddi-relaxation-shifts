"""decomp 1 단계 — 덩어리 하나의 합성 신호 · 확장 좌표 (SPEC_v1 §3 · §1-2). `envs/reference` 에서 돈다.

사용:  python synth.py <exp> <prot> <chunk_id> <outdir>
  - 칸 목록(config.cells)을 프로토콜로 거르고 반복을 펼친 행을 CHUNK 개씩 자른 덩어리 <chunk_id> 를 만든다.
  - outdir 에 data.nii · nodif_brain_mask.nii(압축 없음 — 쓰기 시간) · bvals · bvecs · DWI.scheme(적합기 입력)과
    rows.npz(행 메타 · 겉보기 분율 · 확장 좌표)를 쓴다. 적합값은 읽지 않는다.
구획 신호는 AMICO 2.1.1 합성 클래스의 섬유 방향별 함수(`prereg/scripts/b10_synth.py` 와 같은 함수 · 같은 인자 단위).
"""
from __future__ import annotations

import math
import os
import sys

import nibabel as nib
import numpy as np
from amico.scheme import Scheme
from amico.synthesis import NODDIExtraCellular, NODDIIntraCellular, NODDIIsotropic

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C  # noqa: E402

B0_MAX = 50.0


class Proto:
    """한 프로토콜의 구획 신호 계산기. b>0 인데 방향 벡터가 0 인 행(b = 5)은 x · y · z 세 방향 평균."""

    def __init__(self, name):
        self.name = name
        self.bvals, self.bvecs = C.protocol(name)
        nrm = np.linalg.norm(self.bvecs, axis=1)
        self.b0 = self.bvals <= 0
        self.zero = (self.bvals > 0) & (nrm < 1e-6)
        self.schemes = []
        for e in ((1.0, 0, 0), (0, 1.0, 0), (0, 0, 1.0)) if self.zero.any() else ((1.0, 0, 0),):
            g = self.bvecs.copy()
            g[self.zero | self.b0] = e
            g = g / np.linalg.norm(g, axis=1, keepdims=True)
            sch = Scheme(np.c_[g, self.bvals], b0_thr=0)
            ic, ec, iso = NODDIIntraCellular(sch), NODDIExtraCellular(sch), NODDIIsotropic(sch)
            P = ic.protocol_hr
            args = (P["grad_dirs"], np.squeeze(P["gradient_strength"]), np.squeeze(P["delta"]), np.squeeze(P["smalldel"]))
            self.schemes.append((ic, ec, iso, args))
        self._ic = {}

    def _avg(self, f):
        outs = [f(*s) for s in self.schemes]
        out = outs[0].copy()
        if len(outs) > 1:
            out[self.zero] = np.mean([o[self.zero] for o in outs], axis=0)
        out[self.b0] = 1.0
        return out

    def E_ic(self, dpar, kappa, fdir):
        key = (dpar, kappa, fdir)
        if key not in self._ic:
            d = C.FIBRE_DIRS[fdir]
            self._ic[key] = self._avg(lambda ic, ec, iso, a: np.asarray(
                ic._synth_meas_watson_SH_cyl_neuman_PGSE(np.array([dpar * 1e-9, 0, kappa]), *a, d), float))
        return self._ic[key]

    def E_ec(self, dpar, dperp, kappa, fdir):
        d = C.FIBRE_DIRS[fdir]
        return self._avg(lambda ic, ec, iso, a: np.asarray(
            ec._synth_meas_watson_hindered_diffusion_PGSE(np.array([dpar * 1e-9, dperp * 1e-9, kappa]), *a, d), float))

    def E_iso(self):
        return self._avg(lambda ic, ec, iso, a: np.asarray(iso.get_signal(C.DISO * 1e-3), float))


def noiseless(pr, cell):
    """칸 하나의 무잡음 신호. 반환: (S/S0 정규화 신호, S0 절대값(TE=0 · S_ref=1), ν′, f′_iso)."""
    nu, odi, fiso = (C.TISSUE_E5 if cell.get("dense") else C.TISSUE)[cell["tissue"]]
    if "eta6" in cell:                                          # 3 단계 — η 를 직접(SPEC_v3 §2)
        _, t2, t1 = C.ETA_E6[cell["eta6"]]
        a = C.weights(cell["TE"], t2, cell["TR"], t1)
    else:
        t2 = C.t2_triplet(C.t2_of(cell["t2set"]))
        if np.isnan(cell["TR"]):
            a = C.weights(cell["TE"], t2)
        else:
            a = C.weights(cell["TE"], t2, cell["TR"], C.t1_triplet(C.T1_SETS[cell["t1set"]]))
    nu_app, f_app = C.apparent(nu, fiso, a)
    kp = C.odi_to_kappa(odi)
    dp = cell["dpar"]
    E_ic = pr.E_ic(dp, kp, cell["fdir"])
    E_ec = pr.E_ec(dp, C.dperp(cell["gen"], dp, nu, nu_app), kp, cell["fdir"])
    S = (1 - fiso) * (nu * a[0] * E_ic + (1 - nu) * a[1] * E_ec) + fiso * a[2] * pr.E_iso()
    S0 = (1 - fiso) * (nu * a[0] + (1 - nu) * a[1]) + fiso * a[2]
    return S / S0, S0, nu_app, f_app


def n2_sigma(snr0):
    """Gong 관례: 순수 WM 기준 화소의 S(b=0, TE=98)/σ = SNR₀ (SPEC_v1 §6 · N2_REF)."""
    r = C.N2_REF
    a = np.exp(-r["TE"] / np.array([r["T2_in"], r["T2_en"], 1000.0]))
    S0 = (1 - r["f_iso"]) * (r["f_in"] * a[0] + (1 - r["f_in"]) * a[1]) + r["f_iso"] * a[2]
    return S0 / snr0


def add_noise(S, noise, snr, rng, S0abs=1.0):
    """S: (R, V) 신호. N1/L* — 자기 b0 정규화 신호에 σ = 1/SNR. N2 — 절대 신호(S·S0abs)에 고정 σ."""
    if noise == "N0":
        return S
    if noise == "N2":
        S = S * S0abs
        sig = n2_sigma(snr)
    else:
        sig = 1.0 / snr
    L = {"N1": 1, "N2": 1, "L4": 4, "L8": 8}[noise]
    re = S + rng.normal(0, sig, S.shape)
    acc = re ** 2 + rng.normal(0, sig, S.shape) ** 2
    for _ in range(L - 1):
        acc += rng.normal(0, sig, S.shape) ** 2 + rng.normal(0, sig, S.shape) ** 2
    return np.sqrt(acc)


def ext_coords(S, bvals, bvecs, low):
    """확장 좌표: MD · FA(b0 + 낮은 셸 · 로그 OLS 뒤 S² 가중 1 회) · s_low · s_2000 (SPEC_v1 §1-2)."""
    b0 = bvals <= B0_MAX
    lo = np.abs(bvals - low) <= 50
    hi = np.abs(bvals - 2000.0) <= 50
    S0 = S[:, b0].mean(1)
    S0 = np.maximum(S0, 1e-12)
    s_lo = np.log(np.maximum(S[:, lo].mean(1), 1e-12) / S0)
    s_hi = np.log(np.maximum(S[:, hi].mean(1), 1e-12) / S0)
    rows = b0 | lo
    g, b = bvecs[rows], bvals[rows]
    X = np.c_[np.ones(rows.sum()), -b * g[:, 0] ** 2, -b * g[:, 1] ** 2, -b * g[:, 2] ** 2,
              -2 * b * g[:, 0] * g[:, 1], -2 * b * g[:, 0] * g[:, 2], -2 * b * g[:, 1] * g[:, 2]]
    Y = np.log(np.maximum(S[:, rows], 1e-12))
    beta = np.linalg.lstsq(X, Y.T, rcond=None)[0].T                       # OLS
    w = np.exp(2 * (beta @ X.T))                                          # S² 가중(OLS 예측)
    XtWX = np.einsum("vr,ri,rj->vij", w, X, X)
    XtWy = np.einsum("vr,ri,vr->vi", w, X, Y)
    beta = np.linalg.solve(XtWX, XtWy[..., None])[..., 0]
    D = np.zeros((len(S), 3, 3))
    D[:, 0, 0], D[:, 1, 1], D[:, 2, 2] = beta[:, 1], beta[:, 2], beta[:, 3]
    D[:, 0, 1] = D[:, 1, 0] = beta[:, 4]
    D[:, 0, 2] = D[:, 2, 0] = beta[:, 5]
    D[:, 1, 2] = D[:, 2, 1] = beta[:, 6]
    ev = np.linalg.eigvalsh(D)
    md = ev.mean(1)
    fa = np.sqrt(1.5 * ((ev - md[:, None]) ** 2).sum(1) / np.maximum((ev ** 2).sum(1), 1e-30))
    return np.c_[md * 1e3, fa, s_lo, s_hi].astype(np.float32)               # MD 는 µm²/ms


def main():
    exp, prot, cid, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
    os.makedirs(out, exist_ok=True)
    cl, bounds = C.chunk_rows(exp, prot)
    a, b = bounds[cid]
    pr = Proto(prot)
    rng = np.random.default_rng([20261020, (C.EXPERIMENTS + ("E5", "E6", "E6n")).index(exp), C.PROTOCOLS.index(prot), cid])
    sig, meta = [], {k: [] for k in ("cell", "rep", "nu_app", "f_app")}
    for ci in range(a, b):
        c = cl[ci]
        S, S0abs, nu_app, f_app = noiseless(pr, c)
        R = c["nrep"]
        Sn = add_noise(np.repeat(S[None, :], R, 0), c["noise"], c["snr"], rng, S0abs)
        sig.append(Sn.astype(np.float32))
        meta["cell"] += [ci] * R
        meta["rep"] += list(range(R))
        meta["nu_app"] += [nu_app] * R
        meta["f_app"] += [f_app] * R
    sig = np.concatenate(sig, 0)
    if not np.isfinite(sig).all():
        sys.exit(f"⛔ 비유한 신호 {np.size(sig) - np.isfinite(sig).sum()} 개")
    ext = ext_coords(sig.astype(float), pr.bvals, pr.bvecs, C.low_shell(prot))
    n = len(sig)
    vol = (sig * 1000.0).reshape(n, 1, 1, -1)
    nib.save(nib.Nifti1Image(vol, np.eye(4)), f"{out}/data.nii")
    nib.save(nib.Nifti1Image(np.ones((n, 1, 1), np.uint8), np.eye(4)), f"{out}/nodif_brain_mask.nii")
    np.savetxt(f"{out}/bvals", pr.bvals[None, :], fmt="%g")
    np.savetxt(f"{out}/bvecs", pr.bvecs.T, fmt="%.6f")
    g = pr.bvecs.copy()
    nrm = np.linalg.norm(g, axis=1)
    g[nrm > 0] /= nrm[nrm > 0, None]
    np.savetxt(f"{out}/DWI.scheme", np.c_[g, pr.bvals], fmt="%.6f", header="VERSION: BVECTOR", comments="")
    cols = {k: np.array([cl[i][k] for i in meta["cell"]]) for k in ("TE", "TR", "t2set", "t1set", "tissue", "fdir", "dpar", "snr")}
    for k in ("acq", "eta6", "side", "case", "ref"):
        if k in cl[a]:
            cols[k] = np.array([cl[i][k] for i in meta["cell"]])
    np.savez(f"{out}/rows.npz", exp=exp, prot=prot, chunk=cid, cell=np.array(meta["cell"]), rep=np.array(meta["rep"]),
             gen=np.array([cl[i]["gen"] for i in meta["cell"]]), noise=np.array([cl[i]["noise"] for i in meta["cell"]]),
             nu_app=np.array(meta["nu_app"]), f_app=np.array(meta["f_app"]), ext=ext, **cols)
    print(f"{exp} {prot} 덩어리 {cid}/{len(bounds)}: 행 {n:,} · 볼륨 {sig.shape[1]}", flush=True)


if __name__ == "__main__":
    main()
