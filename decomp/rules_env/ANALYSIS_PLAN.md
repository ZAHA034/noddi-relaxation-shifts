# rules_env 분석 계획 — 계산 전 고정

이 문서는 측정 함수 환경(`decomp/rules_env/`)으로 계산할 분석의 입력 · 단위 · 규칙 · 주 숫자 · 읽기를 **새 양을 계산하기 전에** 고정한다. 이 문서를 담은 커밋의 시각이 «계산 전» 의 기록이다. 저자 결정(2026-11-14): 기준은 재검사에 묶는다 · 실자료 ROI 값 다시 읽기 · 결과를 원고에 넣는다.

## 0. 등급과 범위

- 모든 분석은 **exploratory** 이다(원고 3.3 의 정의 — 실자료의 x · y 값과 Wong 의 TE 기울기, 그리고 같은 시뮬레이션 적합 결과의 방향은 이미 보았다). 원고의 판정(verdict)은 하나도 바꾸지 않는다.
- 시뮬레이션은 **무잡음만** 다룬다. 잡음 바닥 모의는 앞선 탐색 분석(보충 표 S6 No. 28)에서 결론이 무효였고, 시뮬레이터의 백질 FISO 가 실측과 어긋난다(보충 S-M8.5).
- 이미 끝낸 확인(새 양 아님): 새 로더가 기존 보고값을 다시 냄 — Tong 몫 x 0.931 · ODI 0.934 · y 0.844, Wong 선형 몫 x 0.815 · y 0.901, ON-Harmony ODI 재검사 수준 2/19(`rules_env.regress`). 새 합성 · 적합 경로가 저장된 2 단계 적합값을 차 0 으로 다시 냄(`check_synthfit`).

## 1. 공통 규칙(새 문턱 없음)

- **ORS(오프셋 제거 재검사 몫).** 단위 u 의 차이 d_u 에서 같은 묶음의 나머지 단위 d 의 중앙값(자기 제외)을 빼고, |d_u − o_{−u}| ≤ c·√2·s_u 이면 재검사 수준. 몫이 0.90 이상이면(정수 꼴 `analysis.retest_level` 과 같음) «위치 보정으로 재검사 수준». 단위 10 개 미만이거나 뺀 단위가 10 % 를 넘으면 판정 불가. c · s 는 자료마다 원고 3.4 의 규칙 그대로: Tong — t₀.₉₇₅(df) 와 c10 반복의 단위별 SD, ON-Harmony — 1.96 과 OXF1PRI 반복(2 명 × 6 세션)의 단위별 SD, Wong — 1.96 과 반복 계열(TE62R2 · TE132R2)에서 반복 쌍마다 평균 차를 뺀 SD. 척도를 바꾸면 같은 반복에서 그 척도의 SD.
- **ΔTE*(TE 허용폭).** 오프셋 제거 차이 ~ N(0, 2s²(1 + 1/(n−1)) + ΔTE²τ²) 에서 ORS 기대값이 0.90 이 되는 ΔTE. τ 는 단위 사이 TE 기울기의 SD(잡음 분산을 뺌). 근이 70 ms(자료의 TE 폭)를 넘으면 «≥ 70 ms».
- **S2 읽기.** 적합 결과의 사이트 사이 집단 차 변화 ÷ 식 (7)(T1 항 포함)이 모든 해당 칸에서 0.80–1.25 안이면 «식 (7) 이 적합 결과에서 성립». 0.80–1.25 는 평균 생물학적 동등성 실무에서 쓰는 범위를 관례로 가져온 것이다(구간으로 동등성을 판단하는 원리: Schuirmann 1987, J Pharmacokinet Biopharm 15:657–680, doi:10.1007/BF01068419 — 이 문헌이 ±20 % 를 정한 것은 아니다).
- **척도 비교(E4).** 같은 단위에서 두 척도의 ORS 차와 짝지은 부트스트랩 95 % 구간(참여자가 둘 이상이면 참여자를, 아니면 단위를 다시 뽑음 · B 2000 · 시드 0). 그리고 모수 하나짜리 두 모형 «raw 이동 상수» 대 «logit 이동 상수»(d = expit(logit ν_A + c) − ν_A)의 잔차 제곱합.
- **기술(판정 없음).** 비례 편향(Bland & Altman 1999: 차이를 평균에 회귀한 기울기) · 스피어먼 ρ · 순위 뒤집힘 몫(A 에서 c·√2·s 보다 떨어진 단위 쌍 가운데).

## 2. 시뮬레이션(무잡음 · F2 = AMICO 2.1.1 · d∥ 1.7 · d_iso 3.0 µm²/ms)

**S1(a) 저장된 3 단계 적합 결과.** `analyze3.load_table("E6")`(섬유 방향 넷 평균 좌표). 먼저 X9 의 원뿔 안 몫 12 개(생성기 셋 × Q1 – Q4, `decomp/results/stage3/x9.json`)를 넷째 자리까지 다시 낸다 — 어긋나면 멈춤. 그다음 같은 표: 단위 = 기준 조직 210(`config.refs_e5()`), η = L0|iii 와 L1|cv0.03|rho0.0|iii 16 점(17 설정), 쌍 = A 대 Q1(TE 만) · Q4(TR 만) · Q3(방향 구성만) · Q2(모두 다름), 생성기 셋. d_u = x_q − x_A. s = ON-Harmony OXF1PRI 반복으로 낸 x 의 단위별 SD 의 중앙값(27 단위) · c 1.96. 산출: 설정마다 ORS · 비례 편향 · ρ · 순위 뒤집힘. **주 숫자:** 쌍 부류마다 T-in 의 17 설정 가운데 ORS ≥ 0.90 인 몫(다른 생성기는 함께 보고 · T-app 는 명제 2 상 정확한 경우라 양성 대조).

**S1(b) ROI 맞춤 단위 · 새 적합.** Wong: 참여자마다 TE 62 ms 의 27 단위 관측 (NDI, FISO, ODI)를 중심 이완(T2 82 · 51.25 · 1000 ms, T1 1100 · 1100 · 4300 ms, TR 5800 ms)의 식 (3) 역으로 참 (ν, f) 로(ODI 는 그대로 · 참 ν · f 가 (0, 1) 밖이면 그 단위를 빼고 수를 보고) → Wong 카드 8 TE 로 합성 · 적합(생성기 셋 · 섬유 방향 넷 · η 17 설정). 단위마다 x 의 TE 기울기(8 TE OLS)와 8 TE 평균 → τ_fit(단위 사이 SD) · β_fit(기울기와 수준의 공분산 ÷ 수준 분산). ON-Harmony: OXF1PRI 주 분석 피험자의 단위별 중앙값을 OXF1PRI 카드(TE 92 · TR 3600)에서 역변환 → ON-Harmony 카드 여덟(TE · TR 은 기록 범위의 가운데 값 — NOT2ING 의 TR 은 대수 카드 · FISO 부호 규칙과 같은 확정 범위, 설계 S → P2a · V → P2b)으로 합성 · 적합 → 28 쌍의 ORS(s 와 c 는 S1(a) 와 같음).

**S2 명제 4 · 적합 결과.** ROI 조직 = ON-Harmony ROI 맞춤 27 단위. 집단: 대조 · 환자(참 ν −0.03 과 아래 이완). 시나리오: (i) 그림 3 — 대조 T2 82.5 · 52 · 1000, 환자 T2,en 70(T1 없음) (ii) L2 신경돌기 외 — 대조 중심(82 · 51.25 · 1000 · T1 1100 · 1100 · 4300), 환자 T2,en · T1,en × 1.0984 · × 1.0441 (iii) L2 두 구획 — T2,in · T2,en × 1.0984, T1,in · T1,en × 1.0441 (iv) 분율만 — 환자 이완 = 대조. 사람: 사이트 · 집단마다 20 명, 이완은 집단 중심에서 L1 변이(T2 두 구획 CV 0.03 · 상관 0 · 마할라노비스 반지름 2 에서 자름 · 조직 T1 공통 인자 CV 0.04) — 시드 20261115. ROI 값 = 무잡음 적합값(생성기별 · 섬유 방향 넷 평균 좌표) + N(0, s_u²)(s_u = ON-Harmony OXF1PRI 반복의 그 척도 단위별 SD) — 몬테카를로 200 회, 시드 20261116. 사이트: 둘(P2a 체계 · TR 3600 · TE 70 대 100) · 여덟(ON-Harmony 카드). 척도: x(logit) · raw NDI. 방법: 없음 · 대조군 오프셋(대조군 사이트 평균을 단위마다 맞춤) · neuroCombat(위치만 · 위치 + 척도, 집단 공변량) · neuroCombat 을 대조군으로 학습(위치만)해 환자에 적용(`neuroCombatFromTraining` 과 같은 식의 `combat_apply` — 그 함수는 numpy 2 에서 깨진다 · 학습 자료에 학습 추정치를 다시 쓰면 neuroCombat 출력과 같음을 시험). 산출: 사이트별 (환자 − 대조) 차 D_s(단위 평균), 둘 사이트의 변화 D_100 − D_70 와 식 (7) 비율 · 몬테카를로 2.5–97.5 %, 여덟 사이트의 D_s 를 TE 에 회귀한 기울기와 대수 기울기 · 1 자유도 group × TE 항 뒤 남는 사이트 × 집단 분산 · T1 이 다른 시나리오에서 그 잔차와 대수 TR 항의 크기, raw 척도에서 (iv) 의 거짓 사이트 × 집단 변화. **주 숫자:** 둘 사이트 · logit · 시나리오 (i)–(iii) · 위치 방법(없음 · 대조군 오프셋 · ComBat 위치만 · 대조군 학습 위치만) · 생성기 셋의 비율 범위.

**K8 FISO 부호 · 적합 결과.** `fixed/fiso_sign` 의 규칙 그대로(ν 격자 0.30–0.75 · 0.01 간격 · 이완 격자 L4 54 와 L1 48 · 부호는 모든 칸이 같을 때만), Δlogit FISO 를 대수 대신 무잡음 적합값(참 FISO 0.08 · ODI 0.20 · 섬유 방향 넷 평균 · 생성기 셋)으로. 대상: 부호 예측이 있던 8 쌍(장비 쌍 · 카드는 ON-Harmony 카드). 산출: 적합 결과 부호가 대수 부호 · 관측 부호와 같은 쌍 수 — 어느 쪽이든 그대로 보고.

## 3. 실자료(이미 계산한 ROI 값 · 영상 열지 않음)

**E1 Wong(3 명 × 27 단위 × 8 TE).** 참여자마다 단위의 x · y TE 기울기 k̂(8 TE OLS · 반복 계열 제외)와 8 TE 평균 x̄. σ = 평균 차를 뺀 반복 SD(x · y 각각) · var(k̂) 잡음 = σ²/4200 ms² · var(x̄) 잡음 = σ²/8. τ̂ = √max(var(k̂) − σ²/4200, 0), β̂ = cov(k̂, x̄)/(var(x̄) − σ²/8). 합친 값: 참여자 평균을 뺀 k̂ 의 합동 분산에서 잡음 평균을 뺌. ΔTE*(x · y)와 90 % 구간(참여자 묶음 부트스트랩 × 반복 SD 의 √(χ²_df/df) 배 · B 2000 · 시드 0), β̂ 의 95 % 구간(같은 부트스트랩). 모형 점검: 62 대 132 ms 의 ORS 관측값 대 모형 예측. 겉보기 ΔR₂ 의 단위 사이 SD = τ̂_x × 1000 s⁻¹. S1(b) 의 τ_fit · β_fit 와 비교. 민감도: 단위 FISO 수준(8 TE 평균)을 k̂ 에서 회귀로 뺀 뒤의 τ̂ · β̂, 선형 적합 잔차로 낸 σ. **주 숫자:** 합친 ΔTE*_x 와 90 % 구간.

**E2 ON-Harmony(19 판정 가능 장비 쌍).** 척도 여섯(raw NDI · x · raw FISO · logit FISO · y · ODI)마다 쌍 요약(공통 참여자의 B − A 중앙값 · 단위별, `analysis.pair_summary`)의 원래 몫(오프셋 없음)과 ORS, s 는 같은 척도의 OXF1PRI 반복 SD · c 1.96. 단위: 쓰인 세션 어디에서든 로짓이 잘린 단위는 x · y · logit FISO 에서 빼고 수를 보고. 기술: |ΔTE| 와 ORS(x)의 τ_b, E1 의 합친 τ̂ 로 예측한 몫 대 관측 몫.

**E3 Tong(3 명 × 27 단위 × 45 센터 쌍).** 척도 여섯마다 원래 몫(x · y · ODI 는 재현값과 같아야 함)과 ORS(오프셋은 참여자 × 센터 쌍마다 27 단위에서 자기 제외 중앙값) — c = t₀.₉₇₅(df) · s 는 c10 반복의 그 척도 단위별 SD.

**E4 척도 비교.** 주 분석: Wong 62 대 132 ms — x 대 raw NDI, y 대 logit FISO 의 ORS 차와 구간(81 단위 · 참여자 묶음), 그리고 참여자마다 두 모형(raw 이동 상수 · logit 이동 상수)의 잔차 제곱합과 c. 함께: Wong 의 나머지 TE 쌍 27 개, ON-Harmony 19 쌍(척도별 ORS 의 쌍 중앙값 · 쌍마다 x 대 raw NDI 의 차와 단위 부트스트랩 구간), Tong(대조). **주 숫자:** Wong 62–132 의 ORS_x − ORS_rawNDI(구간) · ON-Harmony 19 쌍의 척도별 ORS 중앙값.

## 4. 결과마다 원고에서 바뀌는 것(결과 전에 고정)

| 분석 | 결과 | 원고에 쓰는 것 |
|---|---|---|
| E4 Wong | 구간 > 0 | TE 효과는 ROI 단위에서 logit 척도의 위치 이동이고 raw 척도에서는 아니다 — 식 (6) 과 맞음. 표 2 «What to do» 에 «하모나이즈는 logit NDI 로» |
| E4 Wong | 구간이 0 포함 | 이 ν 범위에서는 자료가 척도를 가르지 못함 — 그대로 씀 |
| E4 Wong | 구간 < 0 | raw 척도가 더 맞음 — 평행이동이 ROI 단위 TE 효과를 설명하지 못함을 5.3 에 씀 |
| E4 ON-Harmony | x · y 의 ORS 중앙값이 raw NDI · logit FISO 이상 | 스캐너 차이도 logit · y 에서 위치 이동에 더 가깝다(기술) — 아니면 그대로 |
| E1 | ΔTE* 와 구간 | 표 2 1 행 근거에 «단위 사이 겉보기 ΔR₂ SD · 위치 보정이 재검사 수준에 머무는 TE 차» · ≥ 70 ms 면 «자료 범위 안에서 깨지지 않음» |
| E1 | β̂ 구간이 0 포함 / 아님 | 수준과 무관(ROI 사이 metric 불변과 맞음) / 이완이 수준과 함께 변함 — 10 ms 당 기울기 변화를 씀 |
| S1 | 쌍 부류별 몫 ≥ 0.90 / 아님 | 그 부류에서 적합 결과의 위치 이동이 재검사 수준 안 / 밖 — 4 절 머리 · 5.3 의 «재지 않음» 문장을 결과로 바꿈 |
| S2 | 비율 모두 0.80–1.25 / 아님 | 식 (7) 이 적합 결과에서 성립 / 벗어난 크기를 씀 — 표 2 6 행 근거 · 2.4 끝 문장 |
| S2 raw | (iv) 의 거짓 변화 크기 | 2.4 의 raw 척도 문장에 적합 결과 크기를 붙임 |
| K8 | 일치 수 | 4.2 의 FISO 부호 실패 옆에 «적합기를 넣은 예측의 부호» 를 씀(설명 여부 그대로) |
| E2 · E3 | 원래 몫 대 ORS | y 의 센터 효과가 위치 이동인지(공격 N3) — 기술로만 |

## 5. 계산 순서와 기록

1. 이 문서와 분석 코드(`analyses.py` · `simrun.py`)를 커밋한다(체인). 2. S1(a) 회귀(X9) → S1 · S2 · K8 → E1 – E4. 3. 결과 json 은 drafts 에 쓰고 drafts git 에 먼저 커밋(장부의 산출 커밋). 4. 원고 v16e · 보충 S-M10 · 표 S17 · 그림 7. 실자료 읽기는 모두 다시 읽기 기록 뒤에 한다(185 건 기록 끝).

## 6. 계산 전 고침(S65 커밋 03631a0 뒤 · 새 양을 계산하기 전 · 2026-11-14)

**고침 1 — S2 의 주 분석.** 2 절에 고정한 설계(사이트 · 집단마다 20 명, 사람마다 L1 이완 변이를 한 번 뽑아 고정하고 몬테카를로로는 ROI 잡음만 바꿈)에서는 사이트 사이 집단 차 변화에 «뽑힌 사람들의 이완 평균 차» 가 섞인다. 대수로 어림한 그 표준편차는 약 0.026 logit(TE 100 ms 에서 사람 한 명의 이동 SD 약 0.068 · 20 명 평균 · 네 칸의 차)으로, 식 (7) 의 값(그림 3 시나리오 −0.148 · L2 두 구획 약 −0.019)과 같은 크기다. 그대로 두면 비율 0.80–1.25 가 적합기와 무관한 표본 오차로 깨진다. 그래서 계산 전에 다음과 같이 고친다.
- **주 분석 = 모집단 모드:** 모든 모의 참여자가 집단 중심 이완(명제 4 의 «모집단 수준 추정» 설정). 몬테카를로 ROI 잡음은 그대로. 비율 읽기(0.80–1.25)와 주 숫자는 이 모드에만.
- **보조 = 원래 설계(L1 변이):** 그대로 계산하되 비교 대상을 같은 사람들의 이완으로 낸 대수 값으로(뽑기 오차가 상쇄됨) · 차이만 보고하고 비율은 내지 않음.
- 코드: `analyses.s2` · `s2_subjects` · `s2_sample_algebra`(시험 둘 더함 — 모집단 모드는 중심 한 사람 · 같은 사람 대수는 이완이 모두 중심이면 식 (7) 과 같음).

## 7. R2 — E4 결과를 본 뒤 정의한 분석(2026-11-15 · 새 양을 계산하기 전 · 모두 exploratory)

E4(Wong 62 대 132 ms 에서 raw NDI 의 오프셋 제거 몫이 x 보다 높음: 차 −0.086, 95 % 구간 −0.148 – −0.037)를 본 뒤 정의했다. 원고의 판정은 바꾸지 않는다. 결과와 상관없이 모두 보고한다.

**공통.** E4 의 규칙 그대로: 참여자마다 27 단위, 단위 차에서 같은 참여자 나머지 단위 차의 중앙값(자기 제외)을 빼고 |·| ≤ 1.96·√2·s(척도마다 그 척도의 반복 SD, `empirical.s_rt_corrected`), 세 참여자를 합친 몫. 몬테카를로는 시드 0. 측정 잡음은 x 에 N(0, σ²)(σ = 그 참여자의 x 반복 SD)로 넣고 raw NDI = expit(x) 로 옮긴다.

**E5 — 단위별 TE 기울기 구조가 E4 를 설명하는가(실자료 ROI 값 다시 읽기 · 영상 열지 않음).**
- 입력: E1 과 같은 양(참여자마다 8 TE 의 x 기울기 k_u · 8 TE 평균 x̄_u · τ̂ · β̂ · σ), 단위의 ν_mid = (raw NDI 62 + 132) / 2.
- 멈춤 문: 같은 로더로 다시 낸 E4 의 x − raw 몫 차가 `e4.json` 의 값과 1e−9 안에서 같지 않으면 멈춘다.
- 1 차 근사(참여자마다): logit 이동의 상대 변화율 R = β̂ / k̄(x 한 단위당), raw 이동이 단위 사이에서 고르려면 R ≈ 2ν̄ − 1(ν̄ = ν_mid 의 평균). raw 이동의 상대 변화율은 R − (2ν̄ − 1). 절댓값이 작은 쪽 척도가 더 고르다고 읽는다.
- 기대 몫(주 숫자): k_u = k̄ + β̂·(x̄_u − 평균) + e_u, e_u ~ N(0, τ_e²), τ_e² = max(τ̂² − β̂²·(var x̄ − σ²/8), 0); 참 x(62) = x̄_u − 35·k_u, x(132) = x̄_u + 35·k_u; 위 공통 규칙으로 ORS_x − ORS_raw 의 평균과 2.5 – 97.5 % (M 2000). 보조: β = 0(같은 τ̂ 를 모두 무작위로).
- 읽기: (가) 기대 차가 E4 의 95 % 구간(−0.148 – −0.037) 안 → «E4 는 E1 의 단위별 이완 구조(수준에 따라 커지는 TE 기울기)로 설명된다» · (나) 밖 → «E1 의 기울기 구조만으로는 E4 를 설명하지 못한다 — 잡음 바닥 · 부분 부피 같은 다른 원인이 남는다». β = 0 판은 기술(읽기 없음).

**E6 — 이완을 공유했다면 E4 가 척도를 가릴 힘이 있었는가(양성 대조 · 저장된 S1b 무잡음 적합값 · 새 적합 없음).**
- 입력: `<작업 폴더>/s1b/fits.npz`(sha256 f0a9c634cacc7b3b171adadcb4e2da2dc7f809dd652f9452094a2485567db384) · `jobs.json`(f9e75337c643e025aaaa6d2b0eed2a920cbcdb3162c91cb614d9edd90607475c). Wong ROI 맞춤 단위 81(S1b 와 같음 · 일 목록이 다시 만든 ROI 맞춤 단위와 다르면 멈춤), TE 62 · 132 ms 의 섬유 방향 평균 x. 이완 설정 하나를 모든 단위가 공유(η 17 × 생성기 3 = 51 칸).
- 잡음과 문턱: 참여자마다 x 의 반복 SD 를 잡음으로, 척도마다 관측 반복 SD 를 문턱으로(E4 와 같음). M 2000(주 칸 T-in · L0|iii) · 500(나머지 50 칸).
- 읽기(주 칸): (가) ORS_x ≥ 0.90 이고 ORS_raw < 0.90 → «E4 는 척도를 가릴 힘이 있었다» · (나) 둘 다 ≥ 0.90 또는 둘 다 < 0.90 → «E4 는 척도에 대해 알려 주는 것이 적다» · (다) ORS_raw ≥ 0.90 이고 ORS_x < 0.90 → «적합기 자체가 raw 를 고르게 만든다 — E4 의 해석이 바뀐다». 나머지 50 칸은 같은 읽기의 칸 수로 기술.

**S3 — raw 척도의 가짜 사이트 × 집단 항(대수 · 자료 없음).**
- 모형: 표 S15 와 같음(`config.weights` TE 항만 · `apparent` 자유수 0 · T2,in 82.5 · T2,en 52 ms · 두 집단이 이완을 공유). 사이트 B TE 70 ms, 사이트 A TE 80 · 90 · 100 ms. 대조 참 ν 0.46(표 S15) · 0.60, 집단 차 Δν −0.03(S2 · 그림 3) · −0.10 · −0.20(큰 집단 차의 예시 · 문헌 값 아님).
- 숫자: raw 집단 차의 사이트 사이 변화 D_A − D_B, logit 의 같은 변화(0 이어야 함 · 대조), 중간 TE 에서 잰 기울기의 근사 Δc·[ν_p(1 − ν_p) − ν_c(1 − ν_c)]. 읽기 없음(크기 기술).
- 기록: 단위 시험을 쓰는 동안(2026-11-15 03:27 KST 무렵 · 이 문서를 커밋하기 전) 실패한 단언의 출력에 한 칸의 raw 값(0.0025 · 대조 0.46 · Δν −0.10 · TE 100 ms 로 보임)이 찍혔다. 그 뒤 시험은 값을 찍지 않게 고쳤다. S3 은 자료가 없는 대수이고 읽기가 없으므로 바뀌는 것은 없다.

**결과마다 원고 문장(결과 전에 고정 · {…} 에 결과 json 의 값).** 놓을 자리: 보충 S-M10(새 소절) · 표 S17 새 부분, 본문 4.3 절에 E5 · E6 한 구절, S3 은 5.4 절 권고의 근거.
- E5 (가): The level dependence of the TE slopes reproduced the scale comparison: with the estimated slopes, the expected difference in offset-subtracted share between x and raw NDI at 62 against 132 ms was {d} ({lo} to {hi} over simulated replicates), within the observed interval, against {d0} without the level dependence (exploratory).
- E5 (나): The level dependence of the TE slopes did not reproduce the scale comparison: the expected difference was {d} ({lo} to {hi}), outside the observed interval, so another cause, such as the noise floor or partial volume, remains (exploratory).
- E5 1 차(공통): To first order, the logit shift grew by {R} of its size per unit of x, against {T} needed for a uniform raw shift ({R − T} on the raw scale).
- E6 (가): Had all ROI units shared relaxation, the same comparison would have separated the scales: at fitter output the expected share was {ox} for x and {or} for raw NDI ({k} of 51 settings with the same reading; exploratory).
- E6 (나): Had all ROI units shared relaxation, the expected shares would have been {ox} for x and {or} for raw NDI, both on the same side of 0.90, so this comparison says little about the scale (exploratory).
- E6 (다): Even with relaxation shared by all ROI units, the fitter output favored raw NDI ({or} against {ox}), which changes the reading of the scale comparison (exploratory).
- S3: On the raw scale, two groups that share relaxation and differ only in true NDI acquire a site-by-group change between TE 70 and 100 ms of {r10} for a difference of 0.10 and {r20} for 0.20 (control NDI 0.46; {r10b} and {r20b} at 0.60), and none on the logit scale (exploratory).

**계산 순서와 기록.** 이 절 · 코드 · 시험을 개발 저장소에 커밋하고 push 한 뒤 `python -B -m rules_env.analyses e5|e6|s3 <출력> <작업>` 를 차례로 돌린다. 결과 json 은 화면에 띄우기 전에 작업 폴더 버전 기록에 커밋한다. 원고에 들어가는 숫자는 장부 행으로 묶는다.

## 8. 구현 오류 고침 — E1 합친 ΔTE*(2026-11-15 · 코드 검토 F2 · E1 결과를 본 뒤 · 다시 계산하기 전)

- **오류.** §1 의 식 N(0, 2s²(1 + 1/(n − 1)) + ΔTE²τ²) 에서 n 은 오프셋(같은 묶음 나머지 단위의 중앙값)을 내는 묶음의 단위 수다. E1 의 묶음은 참여자이고 n = 27 이다. 합친(pooled) 계산은 세 참여자의 81 단위를 한 배열로 넘겨 1 + 1/80 을 썼다(`analyses.e1` · `metrics._te_share`). 참여자별 값과 ON-Harmony 값(묶음 하나 · 27)은 맞다.
- **고친 주 값.** 같은 식 · 같은 부트스트랩(참여자 묶음 · √(χ²/df) · B 2000 · 시드 0)에서 단위마다 그 참여자의 단위 수로 n 을 준다(`analyses.e1_fix` · 결과 `e1_fix.json`). 실행된 값(`e1.json`)은 지우지 않고, 원고 보충 표 S1 (b) 에 이탈로 적는다.
- **멈춤 문.** `e1_fix` 는 참여자별 ΔTE* 를 같은 함수로 다시 내어 `e1.json` 과 1e−9 안에서 같고, 고치지 않은 합친 값(as_run)이 `e1.json` 의 합친 값 · 구간과 같아야 계속한다.
- **민감도(읽기 없음 · 기술).** (i) 합친 τ̂ 의 분산을 N − P 로 나눔(코드 검토 F3) (ii) 오프셋이 평균이 아니라 중앙값이라는 점을 반영한 인자 1 + π/(2(n − 1))(F4) (iii) 둘 다.
- **읽기.** §1 · §4 그대로(ΔTE* < 70 ms 이면 그 값과 90 % 구간). 원고에서는 «about {정수} ms» 와 구간이 바뀌고, 민감도는 보충 S-M10 에 적는다.
- **함께 고친 것(결과 변화 없음).** `combat_apply` 는 공변량 없이 학습한 추정치에만 쓰므로 공변량 추정치를 받으면 멈춘다. `e6` 은 입력 sha256(§7)과 일 목록의 모든 칸을 대조한다(실행된 판은 다섯 칸만 · 사후 대조에서 같음). 변이 시험이 찾은 공백에 시험을 더했다.
