"""Extracts JSON from the four verdict result files in the value area (V4 C3, ON-Harmony pair test, Tong, Wongkornchaovalit) and draws Markdown from it."""
from __future__ import annotations

import hashlib
import json
import os
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rconfig as RC  # noqa: E402

REPO = RC.REPO
WORK = RC.WORK
OUT = REPO / "decomp" / "results" / "remeasure"
OUT_JSON = OUT / "judgements.json"
OUT_MD = OUT / "JUDGEMENTS.md"
ONH_LOG = "prereg/artifacts/decomp/onharmony_image_access_log.json"
EXT_LOG = "prereg/artifacts/decomp/external_data_access_log.json"
INPUTS = {                                  # Name → (path inside the value area, canonical access record)
    "v4_c3": ("onharmony/v4/v4_c3.json", ONH_LOG),
    "onh_pairs": ("onharmony/pairs/onh_pairs.json", ONH_LOG),
    "tong_judge": ("tong/judge/tong_judge.json", EXT_LOG),
    "wong_judge": ("wong/judge/wong_judge.json", EXT_LOG),
}
READ_V4_PASS = "이 설계군 · 이 유형에서 모형에 넣지 않은 스캐너 효과가 판정을 흔들지 않았다"   # Fixed analysis document §2-5 (line 291)
READ_V4_FAIL = "이 유형은 실자료에서 판정 불가"
READ_V4_NODATA = "자료 없음 — 덮음률을 내지 않고 어떤 분모에도 넣지 않는다(a11)"
ACQ = {"match": "matched", "te": "TE only", "tr": "TR only", "tetr": "TE and TR", "tetr_design": "TE, TR and scheme"}
V4_CELL_KEYS = ("acq_type", "blind", "n_subjects", "n_units", "n_units_used", "n_units_dropped_k", "n_excluded_misspec",
                "exceed_rate", "exceed_p95", "gap_detected", "coverage", "coverage_min", "coverage_max", "threshold",
                "threshold_fallback", "passed", "verdict", "snr_neff_dependent", "assumption_dependent",
                "coverage_incl_excluded", "sensitivity_passed", "by_coil_combo", "decomposition")
ALPHA = 0.05                                # §3 test rule 1
POWER_NOTE = "검정력 낮음(19 쌍 — FISO 8 쌍 · 장비 7)"


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def _entries(node):
    if isinstance(node, list) and node and isinstance(node[0], dict) and "utc" in node[0]:
        yield from node
    elif isinstance(node, dict):
        for v in node.values():
            yield from _entries(v)


def logged(rel_file: str, log_path: pathlib.Path) -> tuple[str, str]:
    """(utc, sha256) of the last entry in the access record whose rel_file (relative path from the top of the repository) is given."""
    d = json.loads(log_path.read_text(encoding="utf-8"))
    hits = [(e["utc"], f["sha256"]) for e in _entries(d) for f in e.get("파일", []) if f.get("file") == rel_file]
    if not hits:
        raise SystemExit(f"⛔ 열람 기록 {log_path.name} 에 {rel_file} 항목이 없다")
    return max(hits)


def _v4(v: dict) -> dict:
    def cells(block):
        out = {}
        for key, c in block["cells"].items():
            row = {k: c.get(k) for k in V4_CELL_KEYS if k in c}
            row["reading"] = {"자료 없음": READ_V4_NODATA, "통과": READ_V4_PASS,
                              "실패": READ_V4_FAIL.replace("이 유형은", f"이 유형({ACQ[c['acq_type']]})은")}.get(c["verdict"], c["verdict"])
            out[key] = row
        return out
    m = v["main"]
    return {
        "arm": v["arm"], "seed_numpy": v["seed_numpy"],
        "inputs": {k: v["inputs"][k] for k in ("attr_table_sha256", "test_cover_sha256")},
        "summary": m["summary"],
        "info_none_counts": {comp: len(x) for comp, x in m["info_none"].items()},
        "cells": cells(m),
        "device_concentration_descriptive": {comp: {k: d[k] for k in ("stat", "p", "n_perm")} for comp, d in m["device_concentration_descriptive"].items()},
        "sensitivity": {name: {"summary": v[name]["summary"]} for name in ("with_16794_NOT4GEP001", "without_snr_out_of_train")},
    }


def _sign_read(p: float, what: str) -> str:
    return (f"영가설 기각(유의수준 {ALPHA} 한쪽) — {what}과 맞음 · 비맹검 점검" if p < ALPHA
            else f"영가설 기각되지 않음(유의수준 {ALPHA} 한쪽) — {what}은 지지되지 않음 · 비맹검 점검")


def _onh(d: dict) -> dict:
    m = d["main"]
    oa, ob = m["odi_a"], m["odi_b"]
    if oa["delta_tau"] > 0.5:
        read_a = f"판정 문구 없음 — δ_τ {oa['delta_tau']} > 0.5 → 기술만(검정력 부족)"
    else:
        inside = -oa["delta_tau"] < oa["ci90"][0] and oa["ci90"][1] < oa["delta_tau"]
        read_a = "동등" if inside else "동등이 보이지 않음"
    retest = lambda pp: sorted(k for k, x in pp.items() if x["retest_level"])
    sens = {}
    for name, s in d["sensitivity"].items():
        sens[name] = {"fiso_sign": {k: s["fiso_sign"][k] for k in ("S", "p", "agree", "n_pairs")},
                      "fiso_tau_L1": {k: s["fiso_tau_L1"][k] for k in ("tau", "p", "n_pairs")},
                      "fiso_tau_L4": {k: s["fiso_tau_L4"][k] for k in ("tau", "p", "n_pairs")},
                      "odi_a": {k: s["odi_a"][k] for k in ("tau", "ci90", "n_pairs")},
                      "odi_b": {k: s["odi_b"][k] for k in ("n_pairs", "n_retest_level")}}
    return {
        "rule": d["rule"], "role": m["role"], "power": POWER_NOTE,
        "inputs": d["inputs"],
        "fiso_sign": {**m["fiso_sign"], "reading": _sign_read(m["fiso_sign"]["p"], "부호 예측")},
        "fiso_tau_L1": {**m["fiso_tau_L1"], "reading": _sign_read(m["fiso_tau_L1"]["p"], "크기 순서 예측") + " · p 는 근사"},
        "fiso_tau_L4": {**m["fiso_tau_L4"], "reading": "L1 의 민감도 · p 는 근사"},
        "fiso_null_size": m["fiso_null_size"],
        "odi_a": {**oa, "reading": read_a + " · p 는 근사"},
        "odi_b": {"n_pairs": ob["n_pairs"], "n_retest_level": ob["n_retest_level"], "retest_pairs": retest(ob["per_pair"]),
                  "per_pair": ob["per_pair"], "broken": ob["broken"],
                  "reading": f"{ob['n_pairs']} 쌍 중 {ob['n_retest_level']} 쌍 재검사 수준 — 판정 하나로 묶지 않음(기술 · 깨진 쌍은 획득 차이만)"},
        "odi_b_sensitivity": {name: {"n_pairs": s["n_pairs"], "n_retest_level": s["n_retest_level"], "retest_pairs": retest(s["per_pair"])}
                              for name, s in m["odi_b_sensitivity"].items()},
        "moved_to_not_judgeable": m["moved_to_not_judgeable"],
        "pairs": m["pairs"],
        "sensitivity": sens,
        "sessions": d["sessions"],
        "s_ref_F2": d["s_ref_F2"],
        "clipped": d["clipped"],
    }


def _tong(t: dict) -> dict:
    def block(b):
        ta = lambda x: {k: x[k] for k in ("n_pairs", "n_determined", "n_undetermined", "margin", "undetermined_limit", "p_random", "judgeable", "verdict", "ratio", "p")}
        return {"centers": b["centers"], "s_xy_unit_median": b["s_xy_unit_median"], "df_by_coord": b["df_by_coord"],
                "center_effect_inputs": b["center_effect_inputs"], "T_A": ta(b["T_A"]),
                "T_A_standardized_descriptive": ta(b["T_A_standardized_descriptive"]),
                "T_B": {c: {k: b["T_B"][c][k] for k in ("n_total", "n_den", "n_retest", "ratio", "units_dropped", "judgeable", "verdict")} for c in ("x", "y", "ODI")},
                "T_C": b["T_C"]}
    return {"main": block(t["main"]), "without_c09": block(t["without_c09"]),
            "all_shells_descriptive": block(t["all_shells_descriptive"]), "missing_sessions": t["missing_sessions"]}


def _wong(w: dict) -> dict:
    return {"overall": w["main"]["overall"], "by_subject": w["main"]["by_subject"],
            "without_own_registration": w["without_own_registration"], "own_registration": w["own_registration"],
            "missing_subjects": w["missing_subjects"]}


def extract() -> dict:
    ins, data = {}, {}
    for name, (sub, log) in INPUTS.items():
        p = WORK / sub
        rel = os.path.relpath(p, REPO)
        utc, want = logged(rel, REPO / log)
        got = sha256(p)
        if got != want:
            raise SystemExit(f"⛔ {name}: sha256 이 열람 기록과 다르다(기록 {want[:16]}… · 지금 {got[:16]}…)")
        ins[name] = {"file": rel, "sha256": got, "열람 기록": log, "utc": utc}
        data[name] = json.loads(p.read_text(encoding="utf-8"))
    return {
        "무엇": "⑪ 판정 결과 요약 — 결과 파일 넷에서 스크립트로 옮김(손으로 옮긴 수치 없음) · 판정은 고정 문서 규칙 그대로",
        "생성": "decomp/remeasure/summarize_judgements.py",
        "규칙": "decomp/remeasure/MEASUREMENT_FIXED.md(5f95b18) §2-1 · §2-2 · §2-3 · §2-5 · §2-6 · §3",
        "입력": ins,
        "V4": _v4(data["v4_c3"]),
        "ON-Harmony 쌍 검정": _onh(data["onh_pairs"]),
        "Tong": _tong(data["tong_judge"]),
        "Wongkornchaovalit": _wong(data["wong_judge"]),
    }


def _f(x, n=3):
    return "—" if x is None else (f"{x:.{n}f}" if isinstance(x, float) else str(x))


def render(s: dict) -> str:
    L = []
    w = L.append
    w("# ⑪ 판정 결과 요약")
    w("")
    w(f"{s['무엇']}. 생성: `{s['생성']}`(JSON → 이 MD). 규칙: {s['규칙']}.")
    w("")
    w("| 입력 | 파일(저장소 맨 위 기준) | sha256 | 열람 기록 · utc |")
    w("|---|---|---|---|")
    for k, e in s["입력"].items():
        w(f"| {k} | `{e['file']}` | `{e['sha256']}` | `{e['열람 기록']}` · {e['utc']} |")
    v = s["V4"]
    w("")
    w(f"## V4 — 팔 {v['arm'].upper()}(주 판정 · R11)")
    w("")
    w(f"입력 기록: 속성 표 sha256 `{v['inputs']['attr_table_sha256']}` · 시험 자료 덮음 sha256 `{v['inputs']['test_cover_sha256']}` · numpy 시드 {v['seed_numpy']}. "
      "C7(보조 보고)은 이 요약에 없다 — 결과 파일이 없다(handoff «판정» 의 이탈).")
    w("")
    w("| 성분 | 역할 | «식별됨» 칸 | 덮음률을 낸 칸 | 통과 | 실패 | 판정 불가 | 자료 없음 | 정보 없음(식별 아님) |")
    w("|---|---|---|---|---|---|---|---|---|")
    for comp, x in v["summary"].items():
        w(f"| {comp} | {x['blind']} | {x['tested']} | {x['with_data']} | {x['passed']} | {x['failed']} | {x['not_judgeable']} | {x['no_data']} | {x['info_none']} |")
    w("")
    w("칸마다 따로 판정한다(칸을 모아 비율로 내지 않는다). 미리 정한 읽기 — 통과: «" + READ_V4_PASS + "» · 실패: «" + READ_V4_FAIL + "».")
    w("")
    w("| 칸 | 유형 | 역할 | 피험자 | 단위(쓴 · 전체) | 오지정으로 뺌 | 덮음률 중앙값 [k 최소–최대] | 칸 문턱(대체) | 간극 검출 | SNR/N_eff 의존 | 판정 | 읽기 |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for key, c in v["cells"].items():
        kk = key.replace("|", "\\|")                 # So that «|» in a cell name does not break the MD table
        if c["verdict"] == "자료 없음":
            w(f"| {kk} | {ACQ[c['acq_type']]} | {c['blind']} | 0 | — | — | — | — | — | — | 자료 없음 | {c['reading']} |")
        else:
            w(f"| {kk} | {ACQ[c['acq_type']]} | {c['blind']} | {c['n_subjects']} | {c['n_units_used']} · {c['n_units']} | {c['n_excluded_misspec']} | "
              f"{_f(c['coverage'])} [{_f(c['coverage_min'])}–{_f(c['coverage_max'])}] | {_f(c['threshold'], 2)}{' (대체)' if c['threshold_fallback'] else ''} | "
              f"{'예' if c['gap_detected'] else '아니오'} | {'예' if c['snr_neff_dependent'] else '아니오'} | {c['verdict']} | {c['reading']} |")
    w("")
    w("판정을 바꾸지 않는 보고:")
    for name, x in v["sensitivity"].items():
        parts = [f"{comp} 덮음률을 낸 칸 {y['with_data']} 중 통과 {y['passed']} · 실패 {y['failed']} · 판정 불가 {y['not_judgeable']} · 자료 없음 {y['no_data']}" for comp, y in x["summary"].items() if y["tested"]]
        w(f"- 민감도 {name}: " + " / ".join(parts))
    for comp, x in v["device_concentration_descriptive"].items():
        w(f"- 장비 집중(기술 · 판정에 쓰지 않음) {comp}: 통계량 {x['stat']} · p {_f(x['p'])} · 순열 {x['n_perm']}")
    o = s["ON-Harmony 쌍 검정"]
    w("")
    w("## ON-Harmony 쌍 검정")
    w("")
    w(f"역할: {o['role']} · {o['power']} · 결과 파일 규칙: {o['rule']}.")
    w("")
    w("| 검정 | 통계량 | 장비 순열 p | 효과 크기 · 구간 | 읽기 |")
    w("|---|---|---|---|---|")
    fs, t1, t4, oa, ob = o["fiso_sign"], o["fiso_tau_L1"], o["fiso_tau_L4"], o["odi_a"], o["odi_b"]
    w(f"| FISO 부호 | S = {_f(fs['S'], 0)} · 쌍 {fs['n_pairs']} | 한쪽 {_f(fs['p'])} ({fs['n_perm']}) | 일치 {fs['agree']} / {fs['n_pairs']} = {_f(fs['agree_rate'])} | {fs['reading']} |")
    w(f"| FISO 크기 순서 τ_b — L1 | {_f(t1['tau'])} · 쌍 {t1['n_pairs']} | 한쪽 {_f(t1['p'])} ({t1['n_perm']}) | — | {t1['reading']} |")
    w(f"| FISO 크기 순서 τ_b — L4 | {_f(t4['tau'])} · 쌍 {t4['n_pairs']} | 한쪽 {_f(t4['p'])} ({t4['n_perm']}) | — | {t4['reading']} |")
    w(f"| ODI (a) τ_b | {_f(oa['tau'])} · 쌍 {oa['n_pairs']} | τ = 0 양쪽 {_f(oa['p_tau0_two_sided'])} ({oa['n_perm']}) | 90 % 구간 [{_f(oa['ci90'][0])}, {_f(oa['ci90'][1])}] | {oa['reading']} |")
    w(f"| ODI (b) 재검사 수준 | {ob['n_retest_level']} / {ob['n_pairs']} 쌍 | — | 재검사 수준 쌍: {', '.join(ob['retest_pairs']) or '없음'} | {ob['reading']} |")
    w("")
    w(f"FISO 영 크기(모의 · 고정 문서 §3): {o['fiso_null_size']}.")
    w("")
    w("| 장비 쌍 | 공통 피험자 FISO · ODI | 판정 가능 | Δ logit f_iso | ODI 차 | (b) 재검사 단위 / 값 있는 단위 (필요) | (b) 재검사 수준 |")
    w("|---|---|---|---|---|---|---|")
    for p in o["pairs"]:
        pp = ob["per_pair"].get(p["pair"])
        if pp:
            w(f"| {p['pair']} | {p['n_subjects_fiso']} · {p['n_subjects_odi']} | 예 | {p['fiso']:+.4f} | {p['odi']:+.5f} | {pp['n_retest']} / {pp['n_units']} ({pp['need']}) | {'예' if pp['retest_level'] else '아니오'} |")
        else:
            w(f"| {p['pair']} | {p['n_subjects_fiso']} · {p['n_subjects_odi']} | 아니오(공통 피험자 < 3) | — | — | — | — |")
    w("")
    w("깨진 쌍의 획득 차이(결과 파일 값 그대로 · 방향은 고정 문서 장비 표의 «A~B»):")
    for b in ob["broken"]:
        w(f"- {b['pair']}: ΔTE {b['dTE']:+g} ms · ΔTR {b['dTR_ms']:+g} ms")
    w("")
    w("| ODI (b) 민감도 | 재검사 수준 쌍 | 이름 |")
    w("|---|---|---|")
    for name, x in o["odi_b_sensitivity"].items():
        w(f"| {name} | {x['n_retest_level']} / {x['n_pairs']} | {', '.join(x['retest_pairs']) or '없음'} |")
    w("")
    w("| 민감도 | FISO 부호 S · p · 일치 | FISO τ L1 · p | FISO τ L4 · p | ODI (a) τ_b [90 %] | ODI (b) |")
    w("|---|---|---|---|---|---|")
    for name, x in o["sensitivity"].items():
        a, b1, b4, c, e = x["fiso_sign"], x["fiso_tau_L1"], x["fiso_tau_L4"], x["odi_a"], x["odi_b"]
        w(f"| {name} | {_f(a['S'], 0)} · {_f(a['p'])} · {a['agree']} / {a['n_pairs']} | {_f(b1['tau'])} · {_f(b1['p'])} | {_f(b4['tau'])} · {_f(b4['p'])} | "
          f"{_f(c['tau'])} [{_f(c['ci90'][0])}, {_f(c['ci90'][1])}] | {e['n_retest_level']} / {e['n_pairs']} |")
    ses, sr = o["sessions"], o["s_ref_F2"]
    w("")
    w(f"세션: 뺀 세션 {json.dumps(ses['dropped'], ensure_ascii=False)} · 주 분석에서 뺌 {', '.join(ses['main_excluded'])} · s_ref_F2 세션 {len(ses['s_ref_F2_used'])} · "
      f"자유도 {sr['df']} · 95 % 배수 {_f(sr['multipliers'][0], 4)} – {_f(sr['multipliers'][1], 4)} · 판정 불가로 옮긴 쌍 FISO {len(o['moved_to_not_judgeable']['fiso'])} · ODI {len(o['moved_to_not_judgeable']['odi'])} · "
      f"logit 자름(R25) 세션 {len(o['clipped']['by_session'])} · 쌍 {len(o['clipped']['by_pair'])}.")
    t = s["Tong"]
    w("")
    w("## Tong — 맹검 음성 대조")
    w("")
    w("| 묶음 | T-A 방향 정해진 쌍 / 쌍 · 여유 · 한계 | T-A 판정 | T-B x | T-B y | T-B ODI | T-C |")
    w("|---|---|---|---|---|---|---|")
    for name in ("main", "without_c09", "all_shells_descriptive"):
        b = t[name]
        ta, tb = b["T_A"], b["T_B"]
        cell = lambda c: f"{_f(tb[c]['ratio'])} ({tb[c]['n_retest']} / {tb[c]['n_den']}) {tb[c]['verdict']}"
        w(f"| {name} | {ta['n_determined']} / {ta['n_pairs']} · {_f(ta['margin'])} · {ta['undetermined_limit']} | {ta['verdict']} | {cell('x')} | {cell('y')} | {cell('ODI')} | {b['T_C']} |")
    tm = t["main"]
    w("")
    w(f"주 묶음: 표준화 좌표의 T-A(기술) {tm['T_A_standardized_descriptive']['n_determined']} / {tm['T_A_standardized_descriptive']['n_pairs']} · 여유 {_f(tm['T_A_standardized_descriptive']['margin'])} · "
      f"센터 효과 입력 {json.dumps(tm['center_effect_inputs'], ensure_ascii=False)} · 빠진 세션 {len(t['missing_sessions'])}.")
    g = s["Wongkornchaovalit"]
    w("")
    w("## Wongkornchaovalit — 다중 TE 시연")
    w("")
    w("| 좌표 | 직선 비율 | 단위 수 | 읽기 |")
    w("|---|---|---|---|")
    for c, x in g["overall"].items():
        w(f"| {c} | {_f(x['linear_ratio'])} | {x['n']} | {x['verdict']} |")
    w("")
    w("| 피험자 | x 직선 비율 | y 직선 비율 |")
    w("|---|---|---|")
    for sub, x in g["by_subject"].items():
        w(f"| {sub} | {_f(x['x']['linear_ratio'])} | {_f(x['y']['linear_ratio'])} |")
    w("")
    w(f"따로 정합한 계열: " + " · ".join(f"{sub} {len(x)}" for sub, x in g["own_registration"].items()) + f" · 빠진 피험자 {len(g['missing_subjects'])}.")
    return "\n".join(L) + "\n"


def _write_once(path: pathlib.Path, text: str) -> None:
    if path.exists():
        if path.read_text(encoding="utf-8") != text:
            raise SystemExit(f"⛔ {path.relative_to(REPO)} 이 이미 있고 내용이 다르다 — 덮어쓰지 않는다")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def main(argv: list[str]) -> int:
    if "--render" in argv:
        s = json.loads(OUT_JSON.read_text(encoding="utf-8"))
    else:
        s = extract()
        _write_once(OUT_JSON, json.dumps(s, ensure_ascii=False, indent=1) + "\n")
    _write_once(OUT_MD, render(s))
    for p in (OUT_JSON, OUT_MD):
        print(p.relative_to(REPO), sha256(p))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
