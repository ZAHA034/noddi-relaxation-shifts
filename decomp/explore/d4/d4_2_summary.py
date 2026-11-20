"""Summarizes d4_2_forward.json as Markdown: tables, analytic properties and the ΔTE direction angle, without data."""
import hashlib, json, math, pathlib
O = pathlib.Path(__file__).resolve().parent
J = json.loads((O / "d4_2_forward.json").read_text(encoding="utf-8"))
T2ISO = J["fixed"]["T2ISO"]
ang = lambda t2in, r: math.degrees(math.atan2(1 - t2in / T2ISO, r - 1)) if r != 1 else 90.0
L = ["# D4.2 — 순방향 민감도 표(사후 탐색 · 판정 아님 · 자료 0)", "",
     f"- 스크립트 sha256 {J['script_sha256'][:16]}… · `decomp/stage1/config.py` sha256 {J['config_sha256'][:16]}… · 만든 때 {J['made']}.",
     f"- 고정: T2_iso {T2ISO:.0f} ms · T1_iso {J['fixed']['T1_ISO']:.0f} ms · d_iso {J['fixed']['DISO']} µm²/ms · 기준 조직 ν {J['fixed']['nu']} · f_iso {J['fixed']['f_iso']}.",
     "- NDI · FWF 자리 = 겉보기 분율(ν′ · f′_iso) 과 카드 좌표(x′ = logit ν′ · y′ = logit f′_iso − ln ν′) — AMICO 적합 치우침은 없다.", "",
     "## 1. 카드 좌표는 TE 에 정확히 직선이다(해석 · 수치 확인)", "",
     "- x′ = x − TE·(1/T2_in − 1/T2_en) + T1 항 → 기울기 (r − 1)/T2_in (ms⁻¹) · ν · f · TR 과 무관.",
     "- y′ = y + TE·(1/T2_in − 1/T2_iso) + T1 항 → 기울기 1/T2_in − 1/T2_iso · r · ν · f 와 무관.",
     f"- 수치: x′ 의 직선에서 벗어남 최대 {max(a.get('max_nonlinearity', 0) for a in J['A_x_TE_slope']):.1e} · y′ 은 조직 27 조합 모두 0(62 – 132 ms).",
     "- 그래서 ΔTE 가 카드 점을 옮기는 방향은 각 θ = atan((1 − T2_in/T2_iso)/(r − 1)) 하나로 정해진다(조직과 무관).", "",
     "| T2_in (ms) | r | T2_en (ms) | x′ 기울기 (10 ms 당) | ΔTE 방향 θ (°) |", "|---|---|---|---|---|"]
for a in J["A_x_TE_slope"]:
    rl = f"{a['r']}" if "ref" not in a else f"{a['r']}({a['ref']})"
    L.append(f"| {a['T2_in']} | {rl} | {a['T2_en']} | {a['dx_per_10ms']:+.4f} | {ang(a['T2_in'], a['r'] if 'ref' not in a else a['T2_in'] / a['T2_en']):.1f} |")
th = [ang(t, r) for t in (70.0, 95.0) for r in (1.3, 1.8)]
L += ["", f"- S4 사전(T2_in 70 – 95 · r 1.3 – 1.8)의 θ 범위: {min(th):.1f} – {max(th):.1f}° — 2 단계 REPORT 의 TE 쐐기 48.5 – 72.1° 와 견줄 것(기술).", "",
      "## 2. TR 의 몫(TE 80 ms · r 1.3 · T2_in 82.5) — 1 초당 변화", "", "| T1 비 | T1 공통 c (ms) | TR (ms) | ∂x′/∂TR | ∂y′/∂TR |", "|---|---|---|---|---|"]
for c in J["C_TR"]:
    L.append(f"| {c['T1_ratio']} | {c['T1_c']:.0f} | {c['TR']:.0f} | {c['dx_per_s']:+.5f} | {c['dy_per_s']:+.5f} |")
L += ["", "- T1 비가 1(시나리오 iii)이면 x′ 의 TR 몫은 0 이다. y′ 의 TR 몫은 T1_iso 와 조직 T1 의 차이로 남는다.", "",
      "## 3. 방향 평균 s(b, TE) 의 TE 기울기(10 ms 당) — r 이 기울기를 얼마나 바꾸나", "",
      "| 생성기 | d∥ | r | Wong b700 | Wong b2000 | ON-H b1000 | ON-H b2000 |", "|---|---|---|---|---|---|---|"]
for d in J["D_s_TE"]:
    L.append(f"| {d['generator']} | {d['d_par']} | {d['r']} | {d['Wong_b700_ds_per_10ms']:+.5f} | {d['Wong_b2000_ds_per_10ms']:+.5f} | {d['ON-Harmony_b1000_ds_per_10ms']:+.5f} | {d['ON-Harmony_b2000_ds_per_10ms']:+.5f} |")
L += ["", "- 읽을 것(기술): r 1.0 → 1.8 이 s 의 TE 기울기를 바꾸는 크기를 D3 정답 시험의 잡음 수준(반복 TE62R2 · TE132R2)과 견줘 r 의 식별 가능성을 미리 본다. 판정 아님."]
(O / "d4_2_forward.md").write_text("\n".join(L) + "\n", encoding="utf-8")
print("θ 범위", round(min(th), 1), round(max(th), 1), "| Gong", round(ang(82.0, 82.0 / 52.0), 1), "| E1", round(ang(82.5, 2.64), 1))
print({p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16] for p in sorted(O.iterdir()) if p.is_file()})
