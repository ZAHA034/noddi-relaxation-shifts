"""실측 고정 문서 §3 · §2-1 · §2-2 · §4 — ON-Harmony 쌍 검정(FISO 부호 · FISO τ · ODI (a) · ODI (b)) · s_ref_F2 판정 코드(⑪ 사전 점검 — ⑩ 전에 고정).
값 계산 0 — ⑩ 뒤 ⑪ 이 낸 roi.json 들에 돈다(로더 가드). 규칙은 고정 문서 — 다르면 문서가 정본. 검정 함수는 `decomp/fixed/delta_tau.py` · `analysis.py` 그대로 부른다.

- 세션: judge_v4.session_states 와 같은 규칙(다시도 실패 · 충돌 · 끝 표지 없음은 뺀다 · 일괄 요약이 없으면 멈춤) · 첫 세션 110(§4) ·
  16794_NOT4GEP001 은 주 분석에서 뺀다 · 민감도는 넣는다(R15 · a3).
- 장비 7(NOT2ING 은 한 장비 — 피험자마다 자기 프로토콜의 첫 세션 · a1) · 쌍 21 = delta_tau.PAIRS(장비 이름 순 i < j) · 쌍 값 = D_j − D_i 의
  analysis.pair_summary(단위마다 피험자 중앙값 → 27 단위 중앙값 · R7) · FISO = logit FWF(R2) · ODI = 원척도.
- 판정 가능 쌍: 세션 목록의 공통 피험자 ≥ 3(19 쌍 · §3 표) · 처리 실패를 뺀 뒤 공통 피험자 < 3 이거나 쌍 값이 없으면 판정 불가로 옮긴다(a13) —
  그때만 δ_τ(design_subset) · FISO 영 크기(fiso_sizes)를 바뀐 쌍 집합으로 다시 낸다(설계만 · 같은 시드).
- FISO 부호 = sign_perm(설계 부호 9 쌍 · 관측 부호) · 효과 크기 = 일치 비율(R4) · FISO τ = fiso_tau(최대 |Δ logit f′_iso| L1 · 민감도 L4 · R3) ·
  ODI (a) = τ_b · 90 % 구간(interval90_subset) · τ = 0 의 장비 순열 p(양쪽 — [A 제안 a16]) · δ_τ 0.57 > 0.5 면 기술만(R5) ·
  ODI (b) = 단위마다 |Δ| ≤ 1.96 √2 s · 값과 s 가 있는 단위 중 0.90(retest_level · R6) · s = s_ref_F2(OXF1PRI 반복 2 명 × 6 · analysis.s_ref) ·
  민감도 = 옛 s_ref(5262f43) · s_ref_F2 × 위 끝 배수(sd_multipliers · 남은 세션의 자유도).
- 민감도(같은 검정을 다시): 원척도 ΔFISO(R2) · 16794 넣음(R15) · NOT2ING-98 피험자만(예측도 그 행 · a1) · 위험 ROI 4 종 뺌(R24 · OSF §6-2) ·
  옛 «0 아닌 화소» 표(R13 — roi.json «민감도_옛_-M»).
- 보고: logit 자름 단위 수 — 세션 · 쌍별(R25) · 뺀 세션과 까닭 · 판정 불가로 옮긴 쌍 · 다중 비교 보정 없음(R9).
- 결과 = WORK/onharmony/pairs/onh_pairs.json — 한 번만 쓰고 열람 기록에 sha256(judge_v4 와 같은 방식).

    CUDA_VISIBLE_DEVICES= uv run --frozen python decomp/remeasure/judge_onh.py        # ⑩ · ⑪ 뒤
"""
from __future__ import annotations

import hashlib
import json
import math
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "fixed"))
sys.path.insert(0, str(HERE.parents[0] / "stage1"))
import analysis as AN  # noqa: E402
import delta_tau as DT  # noqa: E402
import fiso_sign as FS  # noqa: E402
import rconfig as R  # noqa: E402

DEVICES, PAIRS = DT.DEVICES, DT.PAIRS
RISK_ROI = ("Fornix (column and body of fornix)", "External capsule", "Cingulum (hippocampus)", "Tapetum")   # R24 · OSF §6-2(§177-6)
DELTA_TAU_FIXED = 0.57                                                    # 고정 문서 §3 — 19 쌍 · 피험자 가중 x · r 격자 최대
FISO_SIZE_FIXED = "모의 4,000 · 시드 20261026: S 0.026 – 0.043 · τ 0.044 – 0.053 / 200,000: S 0.026 – 0.042 · τ 0.046 – 0.052(고정 문서 §3)"
OLD_SREF = HERE.parents[1] / "prereg" / "artifacts" / "unseal" / "s_ref.frozen.json"     # 5262f43 — ODI 축만 쓴다
OLD_SREF_SHA = "39e987aba407a91e6c7b064f052e80e4c610f5c7bda021b51ed4fb33da977849"         # 고정 문서 §2-2 — 읽기 전에 대조(동결 목록 밖 파일)


def pair_name(k: int) -> str:
    i, j = PAIRS[k]
    return f"{DEVICES[i]}~{DEVICES[j]}"


def units_of(rois: dict, table: str) -> list[str]:
    """27 단위 이름 — roi.json 표의 키(roi.py unit_table 과 같은 순서 · 처음 본 세션의 것)."""
    for r in rois.values():
        return list(r[table])
    return []


def device_values(by_unit: dict, rois: dict, feature: str, scale: str, table: str, not2ing: str = "both") -> dict[str, dict[str, dict]]:
    """{장비: {피험자: {단위: 값}}} — by_unit = judge_v4.usable_sessions 꼴({카드 단위: {피험자: 세션}}) · NOT2ING 은 두 프로토콜을 한 장비로
    (not2ing «98» 이면 NOT2ING-98 피험자만 — 민감도)."""
    out = {d: {} for d in DEVICES}
    for unit, subs in by_unit.items():
        if unit == "NOT2ING-89" and not2ing == "98":
            continue
        dev = DT.device_of(unit)
        for sub, ses in subs.items():
            out[dev][sub] = AN.unit_values(rois[f"{sub}_{ses}"], feature, scale, table)
    return out


def pair_values(vals: dict, units) -> list[dict]:
    """PAIRS 순 — 쌍 값 = D_j − D_i(analysis.pair_summary · B − A)."""
    return [AN.pair_summary(vals[DEVICES[i]], vals[DEVICES[j]], units) for i, j in PAIRS]


def keep_pairs(cnt: dict, summ: list[dict]) -> tuple[list[int], list[str]]:
    """판정 가능 쌍 — 세션 목록의 공통 피험자 ≥ 3 이고, 실패를 뺀 뒤에도 공통 피험자 ≥ 3 · 쌍 값이 있을 때(a13 — 값의 크기는 보지 않는다)."""
    keep, moved = [], []
    for k, (i, j) in enumerate(PAIRS):
        if not cnt[frozenset((DEVICES[i], DEVICES[j]))]["judgeable"]:
            continue
        if summ[k]["judgeable"] and summ[k]["value"] is not None:
            keep.append(k)
        else:
            moved.append(pair_name(k))
    return keep, moved


def odi_a(x: np.ndarray, y: np.ndarray, keep: list[int], delta_tau: float) -> dict:
    """ODI (a) — τ_b · 90 % 구간(순열 역산) · τ = 0 의 장비 순열 p(양쪽 [A 제안 a16]) · δ_τ > 0.5 면 기술만 · 아니면 구간이 (−δ_τ, δ_τ) 안이면 «동등»."""
    idx = DT.device_perm_index()
    t, lo, hi = DT.interval90_subset(x, y, keep, idx)
    T = DT.perm_taus_x(x, y, keep, idx)
    desc = delta_tau > 0.5
    return dict(tau=t, ci90=[lo, hi], p_tau0_two_sided=float((np.abs(T) >= abs(t) - 1e-12).mean()), n_perm=len(T), n_pairs=len(keep),
                delta_tau=delta_tau, descriptive_only=desc, equivalent=(None if desc else bool(lo > -delta_tau and hi < delta_tau)))


def odi_b(summ: list[dict], keep: list[int], s: dict, rows_by_pair: dict) -> dict:
    """ODI (b) — 쌍마다 retest_level · «19 쌍 중 몇 쌍» · 깨진 쌍의 획득 차이 칸(ΔTE · ΔTR — 새 설명 0)."""
    per = {pair_name(k): AN.retest_level(summ[k]["by_unit"], s) for k in keep}
    broken = [dict(pair=p, dTE=rows_by_pair[p]["dTE"], dTR_ms=rows_by_pair[p]["dTR_ms"]) for p, r in per.items() if not r["retest_level"]]
    return dict(n_pairs=len(per), n_retest_level=sum(r["retest_level"] for r in per.values()), broken=broken, per_pair=per)


def battery(cnt: dict, rois: dict, by_unit: dict, units, table: str, fiso_scale: str, not2ing: str, keep0: list[int], main_design: bool,
            s_main: dict, s_sens: dict, n_sim_delta: int, n_sim_size: int) -> dict:
    """한 줄기(주 분석 또는 민감도 하나) — 관측 → 판정 가능 쌍 → FISO 부호 · FISO τ(L1 · L4) · ODI (a) · ODI (b).
    main_design = 주 분석과 같은 설계(세션 목록의 공통 피험자 · 예측)인가 — 그렇고 판정 가능 쌍이 keep0 그대로면 고정 문서의 δ_τ · FISO 영 크기를 쓰고,
    아니면(a13 · NOT2ING-98 만) 같은 규칙 · 같은 시드로 다시 낸다(설계만). ODI (b) 의 s 는 모든 줄기에서 s_ref_F2(«주»)와 그 민감도."""
    fv = device_values(by_unit, rois, "FWF", fiso_scale, table, not2ing)
    ov = device_values(by_unit, rois, "ODI", "raw", table, not2ing)
    fs, os_ = pair_values(fv, units), pair_values(ov, units)
    keep_f, moved_f = keep_pairs(cnt, fs)
    keep_o, moved_o = keep_pairs(cnt, os_)
    rows = FS.device_rows(cnt)
    rows_by_pair = {}
    for r in rows:
        a, b = r["pair"].split("~")
        i, j = sorted((DEVICES.index(a), DEVICES.index(b)))
        rows_by_pair[f"{DEVICES[i]}~{DEVICES[j]}"] = r
    s = DT.design_signs(cnt)
    o = np.array([np.sign(v["value"]) if v["value"] is not None else np.nan for v in fs])
    yf = np.array([abs(v["value"]) if v["value"] is not None else np.nan for v in fs])
    yo = np.array([abs(v["value"]) if v["value"] is not None else np.nan for v in os_])
    x_odi = DT.pair_predictions_weighted(cnt, "L1")
    redo_f = not main_design or sorted(keep_f) != sorted(keep0)
    redo_o = not main_design or sorted(keep_o) != sorted(keep0)
    if redo_o:                                                          # a13 — 쌍 집합이 바뀌면 같은 규칙으로 다시(설계만)
        dsub = DT.design_subset(keep_o, x_odi, n_sim=n_sim_delta)
        dt = dsub["delta_tau"] if dsub["delta_tau"] is not None else math.inf
    else:
        dsub, dt = None, DELTA_TAU_FIXED
    x1 = DT.magnitude_vector(rows, "fiso_max_abs_dlogit_L1")
    sizes = DT.fiso_sizes(s, x1, keep_f, n_sim=n_sim_size) if redo_f else FISO_SIZE_FIXED
    return dict(
        pairs=[dict(pair=pair_name(k), fiso=fs[k]["value"], odi=os_[k]["value"], n_subjects_fiso=fs[k]["n_subjects"], n_subjects_odi=os_[k]["n_subjects"],
                    judgeable_fiso=k in keep_f, judgeable_odi=k in keep_o) for k in range(len(PAIRS))],
        moved_to_not_judgeable=dict(fiso=moved_f, odi=moved_o),
        fiso_sign=DT.sign_perm(s, o, keep_f),
        fiso_tau_L1=DT.fiso_tau(x1, yf, s, keep_f),
        fiso_tau_L4=DT.fiso_tau(DT.magnitude_vector(rows, "fiso_max_abs_dlogit_L4"), yf, s, keep_f),
        fiso_null_size=sizes,
        odi_a=odi_a(x_odi, yo, keep_o, dt), odi_a_delta_tau_redone=dsub,
        odi_b=odi_b(os_, keep_o, s_main, rows_by_pair),
        odi_b_sensitivity={k: odi_b(os_, keep_o, v, rows_by_pair) for k, v in s_sens.items()},
        role="FISO 부호 · τ = 비맹검 점검 · ODI (a) · (b) = 기술(§3 역할 표) · 보정 없음(R9)")


def s_ref_f2(work: pathlib.Path, states: dict, rois: dict, units) -> tuple[dict, int, list[str]]:
    """s_ref_F2 — OXF1PRI 반복 2 명 × 6 세션(R.SREF_SUBJECTS · SREF_SESSIONS) 중 쓸 수 있는 세션의 ODI(원척도 · «주») → analysis.s_ref."""
    import judge_v4 as V4  # noqa: E402
    vals, used = {}, []
    for sub in R.SREF_SUBJECTS:
        for ses in R.SREF_SESSIONS:
            key = f"{sub}_{ses}"
            if states.get(key) == V4.USABLE:
                vals.setdefault(sub, []).append(AN.unit_values(rois[key], "ODI", "raw", "주"))
                used.append(key)
    s, df = AN.s_ref(vals, units)
    return s, df, used


def read_rois(work: pathlib.Path, keys) -> dict:
    root = work / "onharmony"
    return {k: json.loads((root / k / "roi.json").read_text(encoding="utf-8")) for k in keys}


def run(work: pathlib.Path = R.WORK, n_sim_delta: int = 1000, n_sim_size: int = 4000) -> tuple[dict, str]:
    """⑩ · ⑪ 뒤 — 주 분석 + 민감도 다섯 → 결과 파일(한 번만) · sha256 · 열람 기록."""
    import guard as G  # noqa: E402
    import judge_v4 as V4  # noqa: E402
    basis = G.require()
    sref_keys = [(a, b) for a in R.SREF_SUBJECTS for b in R.SREF_SESSIONS]
    states = V4.session_states(work, list(dict.fromkeys(R.onh_sessions() + sref_keys)))
    by_unit, dropped = V4.usable_sessions(work)
    main_by_unit = {u: {s: v for s, v in d.items() if (s, v) not in R.MAIN_EXCLUDE_ONH} for u, d in by_unit.items()}
    rois = read_rois(work, [k for k, v in states.items() if v == V4.USABLE])
    units = units_of(rois, "주")
    cnt = AN.common_counts(R.onh_sessions())
    keep0 = [k for k, (i, j) in enumerate(PAIRS) if cnt[frozenset((DEVICES[i], DEVICES[j]))]["judgeable"]]
    s_main, df, used = s_ref_f2(work, states, rois, units)
    lo_m, hi_m = AN.sd_multipliers(df) if df > 0 else (None, None)
    raw = OLD_SREF.read_bytes()
    if hashlib.sha256(raw).hexdigest() != OLD_SREF_SHA:
        raise RuntimeError(f"⛔ {OLD_SREF.name} 의 sha256 이 고정 문서 §2-2 의 값과 다르다 — 옛 s_ref 민감도를 내지 않고 멈춘다")
    old = json.loads(raw.decode("utf-8"))["ROI별 s_ref"]["odi"]
    s_sens = {"옛_s_ref_5262f43": old}
    if hi_m is not None:
        s_sens[f"s_ref_F2_x{hi_m:.4f}"] = {u: (v * hi_m if v is not None else None) for u, v in s_main.items()}
    sess98 = [(a, b) for a, b in R.onh_sessions() if not (R.device_of(b) == "NOT2ING" and a in R.NOT2ING_TE89)]
    cnt98 = AN.common_counts(sess98)
    keep98 = [k for k, (i, j) in enumerate(PAIRS) if cnt98[frozenset((DEVICES[i], DEVICES[j]))]["judgeable"]]
    kw = dict(n_sim_delta=n_sim_delta, n_sim_size=n_sim_size, s_main=s_main, s_sens=s_sens)
    res = dict(
        rule="실측 고정 문서 §3 · §2-1 · §2-2 · §4(⑩ 판) · 코드 judge_onh.py", inputs=dict(old_s_ref_sha256=OLD_SREF_SHA),
        sessions=dict(dropped=dropped, main_excluded=sorted(f"{a}_{b}" for a, b in R.MAIN_EXCLUDE_ONH), s_ref_F2_used=used),
        s_ref_F2=dict(by_unit=s_main, df=df, multipliers=[lo_m, hi_m]),
        clipped=dict(by_session={k: AN.clipped_units(r) for k, r in rois.items() if AN.clipped_units(r)}),
        main=battery(cnt, rois, main_by_unit, units, "주", "logit", "both", keep0, True, **kw),
        sensitivity=dict(
            fiso_raw_scale=battery(cnt, rois, main_by_unit, units, "주", "raw", "both", keep0, True, **kw),
            with_16794_NOT4GEP001=battery(cnt, rois, by_unit, units, "주", "logit", "both", keep0, True, **kw),
            not2ing_98_only=battery(cnt98, rois, main_by_unit, units, "주", "logit", "98", keep98, False, **kw),
            without_risk_roi=battery(cnt, rois, main_by_unit, [u for u in units if u not in RISK_ROI], "주", "logit", "both", keep0, True, **kw),
            old_nonzero_table=battery(cnt, rois, main_by_unit, units_of(rois, "민감도_옛_-M"), "민감도_옛_-M", "logit", "both", keep0, True, **kw)))
    res["clipped"]["by_pair"] = clipped_by_pair(main_by_unit, rois)
    path = work / "onharmony" / "pairs" / "onh_pairs.json"
    sha = V4._write_once(path, res)
    V4._log("ON-Harmony 쌍 검정 — 결과 onh_pairs.json", path, basis)
    return res, sha


def clipped_by_pair(by_unit: dict, rois: dict) -> dict[str, int]:
    """R25 — 쌍마다 공통 피험자의 두 세션 중 하나라도 logit 자름에 든 (피험자 · 단위) 수."""
    sess = {d: {} for d in DEVICES}
    for unit, subs in by_unit.items():
        for sub, ses in subs.items():
            sess[DT.device_of(unit)][sub] = f"{sub}_{ses}"
    out = {}
    for k, (i, j) in enumerate(PAIRS):
        A, B = sess[DEVICES[i]], sess[DEVICES[j]]
        n = sum(len(set(AN.clipped_units(rois[A[s]])) | set(AN.clipped_units(rois[B[s]]))) for s in set(A) & set(B))
        if n:
            out[pair_name(k)] = n
    return out


if __name__ == "__main__":
    _, sha = run()
    print(R.WORK / "onharmony" / "pairs" / "onh_pairs.json", sha)
