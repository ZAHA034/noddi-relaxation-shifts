"""새 실측(⑨-a) — 표준 공간 지도 → 골격 ROI 48 라벨 → 27 단위(옛 등록 규칙: prereg/scripts/lockd/roi_extract.sh:14–28 · measure_high_odi.py:8–9 · OSF EN :247 · :249).

라벨 평균의 화소 집합:
  주(검토자 조건부 받음 · 중계 ZZZZe-2 ②) = 라벨 ∩ 골격(문턱 2000) ∩ 유효 마스크의 모든 화소(0 포함).
       유효 마스크 = 개인 뇌 마스크의 표준 공간 변환(nn) ∩ 적합 유효 지도(F2 유효 · 확장 좌표 유효 — trilinear 로 옮겨 ≥ 0.99).
       뇌 밖 · 정합 밖 · 적합 실패 화소는 값이 0 인지가 아니라 이 마스크로 뺀다. 마스크 안의 진짜 0(예: FWF = 0)은 넣는다 —
       AMICO 는 백질에서 FWF 를 정확히 0 으로 내는 일이 잦아 0 을 빼면 평균이 위로 치우친다. 일곱 지도가 같은 화소 집합을 쓴다(C3 은 C7 의 앞 세 좌표).
  민감도 = 옛 규칙(fslstats -K … -M — 0 이 아닌 화소만 · measure_high_odi.py:8).
27 단위 = 좌우 라벨 값의 단순 평균(가중 없음 · OSF EN :249) · 화소가 없는 라벨은 뺀다(measure_high_odi.py unit_table 과 같다).
"""
from __future__ import annotations

import pathlib
import re
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import rconfig as R  # noqa: E402

N_LABELS = 48


def _labels(xml: pathlib.Path) -> dict[int, str]:
    return {int(m.group(1)): m.group(2).strip() for m in re.finditer(r'<label index="(\d+)"[^>]*>([^<]+)</label>', xml.read_text(errors="replace"))}


def label_units(xml: pathlib.Path = R.ATLAS_XML) -> list[str]:
    """라벨 1..48 → 27 단위 이름(좌우 꼬리표 제거) — prereg/scripts/lockd/_harness.py `_label_units` 와 같은 규칙."""
    lab = _labels(xml)
    return [re.sub(r"\s+[LR]$", "", lab[v]).strip() for v in range(1, N_LABELS + 1)]


def unit_order(xml: pathlib.Path = R.ATLAS_XML) -> list[str]:
    """27 단위의 순서(라벨 번호 순 첫 등장) — `_harness.roi_names()` 와 같다."""
    out = []
    for nm in label_units(xml):
        if nm not in out:
            out.append(nm)
    return out


def skeleton_mask(skel: np.ndarray) -> np.ndarray:
    """fslmaths skel -thr 2000 -bin 과 같다(roi_extract.sh:14)."""
    return skel >= R.SKELETON_THR


EXT = ("MD", "FA", "s_low", "s_2000")
VALID_MIN = 0.99                                                      # 표준 공간 유효 지도(trilinear)가 이보다 작으면 무효 화소가 섞인 것


def valid_mask(brain: np.ndarray | None = None, f2_valid: np.ndarray | None = None, ext_valid: np.ndarray | None = None) -> np.ndarray | None:
    """주 규칙의 유효 마스크(중계 ZZZZe-2 ②) — 개인 뇌 마스크(표준 공간 · nn > 0.5) ∩ F2 유효 ≥ 0.99 ∩ 확장 좌표 유효 ≥ 0.99. 셋 다 없으면 None."""
    parts = [m for m in (None if brain is None else brain > 0.5,
                         None if f2_valid is None else f2_valid >= VALID_MIN,
                         None if ext_valid is None else ext_valid >= VALID_MIN) if m is not None]
    if not parts:
        return None
    out = parts[0].copy()
    for m in parts[1:]:
        out &= m
    return out


def label_means(maps: dict[str, np.ndarray], atlas: np.ndarray, skel: np.ndarray, valid: np.ndarray | None = None,
                nonzero: bool = False, ext_valid: np.ndarray | None = None) -> dict:
    """라벨 1..48 마다 지도 평균과 화소 수.
    주: valid(= valid_mask(...)) 안의 모든 화소 — 0 포함 · 일곱 지도가 같은 화소.
    민감도: nonzero=True 면 그 지도에서 0 이 아닌 화소만(옛 -M 규칙 · 지도마다 따로) · 확장 좌표(MD · FA · s_low · s_2000)는 ext_valid ≥ 0.99 인 화소만."""
    base = skeleton_mask(skel) & (atlas > 0)
    if valid is not None:
        base &= valid
    out = {"n": np.zeros(N_LABELS, int)}
    lab = atlas.astype(int)
    for k, v in maps.items():
        means = np.full(N_LABELS, np.nan)
        cnt = np.zeros(N_LABELS, int)
        for i in range(1, N_LABELS + 1):
            sel = base & (lab == i)
            if ext_valid is not None and k in EXT:
                sel &= ext_valid >= VALID_MIN
            if nonzero:
                sel &= v != 0
            vals = v[sel]
            cnt[i - 1] = vals.size
            if vals.size:
                means[i - 1] = float(vals.mean())
        out[k] = means
        out[f"n_{k}"] = cnt
    out["n"] = np.array([int((base & (lab == i)).sum()) for i in range(1, N_LABELS + 1)])
    return out


def unit_table(lm: dict, features=R.FEATURES, xml: pathlib.Path = R.ATLAS_XML) -> dict[str, dict]:
    """27 단위 = 좌우 라벨 값의 단순 평균 — 값이 없는(nan) 라벨은 뺀다."""
    units = label_units(xml)
    out = {}
    for u in unit_order(xml):
        labs = [i for i, nm in enumerate(units, start=1) if nm == u]
        row = {"라벨": labs, "화소": int(sum(lm["n"][i - 1] for i in labs))}
        for f in features:
            vals = [lm[f][i - 1] for i in labs if np.isfinite(lm[f][i - 1])]
            row[f] = float(np.mean(vals)) if vals else None
        out[u] = row
    return out


def extract(std_maps: dict[str, pathlib.Path], brain_std: pathlib.Path | None, valid_std: pathlib.Path | None = None,
            f2valid_std: pathlib.Path | None = None) -> dict:
    """표준 공간 지도 파일들 → {«주»: 27 단위 표, «민감도_옛_-M»: 27 단위 표, «라벨_주»: 라벨 평균(주)}."""
    import nibabel as nib
    load = lambda p: None if p is None else np.asarray(nib.load(str(p)).dataobj, dtype=float)  # noqa: E731
    atlas = np.asarray(nib.load(str(R.ATLAS)).dataobj)
    skel = np.asarray(nib.load(str(R.FSLDIR / (R.SKELETON + ".nii.gz"))).dataobj)
    maps = {k: load(p) for k, p in std_maps.items()}
    brain, extv, f2v = load(brain_std), load(valid_std), load(f2valid_std)
    main = label_means(maps, atlas, skel, valid_mask(brain, f2v, extv))
    old = label_means(maps, atlas, skel, None, nonzero=True, ext_valid=extv)
    tolist = lambda d: {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in d.items()}  # noqa: E731
    return {"주": unit_table(main, tuple(maps)), "민감도_옛_-M": unit_table(old, tuple(maps)), "라벨_주": tolist(main)}
