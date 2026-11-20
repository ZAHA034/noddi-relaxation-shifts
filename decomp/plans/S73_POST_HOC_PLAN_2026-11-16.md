# S73 방향 메시지: 원고를 더 단단하게 — 문장 고칠 것과 계산 전에 고정하는 분석 넷

> 이 문서는 검토 세션의 방향 메시지 S73(2026-11-16, 저자 확인 뒤 전달)을 받은 그대로 옮긴 것이다. 저장소 규칙에 따라 사람 이름을 «저자» 로 바꾼 것 말고는 고치지 않았다. §2–§5 의 정의 · 판독 · 결과별 문장은 이 커밋으로 고정하며, 이 커밋보다 앞서 §2–§5 의 어떤 양도 계산하지 않았다. 구현에서 이 문서가 열어 둔 선택은 계산 전에 따로 커밋하는 구현 메모(`decomp/plans/S73_IMPLEMENTATION_2026-11-16.md`)에 적는다.

보낸 쪽: 검토 세션(in-project AI reading, 외부 심사 아님). 저자 확인 후 A에게 전달.
대상: v16k 투고판 본문(33쪽, sha256 `1c6456bf…`)과 보충(90쪽, `028b4cc2…`).
표시: [저자 결정] / [검토자 제안] / [저자 결정 필요].
S72의 형식 항목과는 겹치지 않는다. 함께 진행한다.

## 0. 이 메시지의 목적과 원칙

- **목적:** [저자 결정] 원고의 엄밀성과 최종 게재 가능성을 높인다. 아래 §1(문장)과 §2–§5(분석 넷)를 진행한다.
- **원칙:**
  1. **정의 먼저.** §2–§5의 정의, 판독, 결과별 문장을 이 문서 그대로 저장소에 커밋한 뒤 계산한다. 커밋 해시와 시각을 기록한다. 등급은 모두 post hoc이다(판정 뒤에 정의, 계산 전에 고정).
  2. **결과가 나빠도 그대로 보고한다.** 고정 판정(표 1, 보충 표 S5)은 바뀌지 않는다. 새 결과는 고정 결과 옆에 나란히 둔다.
  3. **실패 대조.** 검사마다 일부러 틀린 입력을 넣어, 검사가 그것을 잡는지 먼저 보인다. 대조가 기대대로 움직이지 않으면 그 분석은 멈추고 보고한다.
  4. **CPU만 쓴다.** 새 자료는 내려받지 않는다. 이미 처리한 파일만 쓴다. 새로 내려받거나 새로 읽어야 할 자료가 생기면 멈추고 저자에게 묻는다.
- **실행 순서:** §1은 바로 진행한다. 분석은 §3 → §4 → §2 → §5 순서로 한다. §5의 판독은 §2의 결과를 쓴다.

## 1. 문장 고칠 것 (계산 없음) [검토자 제안]

### 1a. 5.1절의 과장

- 지금 문장:
  > Between regions of one brain, the condition of shared relaxation failed (exploratory).
- 문제: 4.3절은 원인 후보(영역별 구획 T2, 잡음 바닥)를 "분리하지 못했다"고 쓴다. 관찰한 것은 "공통 TE 이동 하나가 영역 사이에 성립하지 않았다"까지다. §2에서 보듯 각 계열을 따로 정합한 것도 후보다.
- 고친 문안(잠정). 최종 문안은 §5.4의 표에 따라 정한다.
  > Between regions of one brain, one common TE shift did not hold (exploratory).
- 같은 단락의 4.3절 문장 «Both results concern a shift shared by regions with different relaxation»은 «…shared by regions that may differ in relaxation»으로 바꾼다.
- 4.3절과 5.4절의 척도 기준 문장은 조건문(«regions with different relaxation»)이므로 그대로 둔다.
- 본문, 보충, README, 커버레터에서 영역 사이 이완 차이를 단정하는 표현을 모두 찾아 목록으로 보고한다(예: «shared relaxation failed», «relaxation differs between regions»). 고치는 것은 그 목록을 본 뒤에 한다.

### 1b. 4.3절의 정합 문장

- 지금 문장:
  > Every series other than the reference failed the alignment check and was registered on its own FA, so the analysis without separately registered series could not be computed
- 문제: 독자는 "자료가 어긋나 있었다"로 읽는다. 그러나 보충 표 S10(d)에 따르면 같은 TE의 반복 계열도 기준 아래였다(r 0.003–0.21). 검사가 정렬된 계열과 어긋난 계열을 구분하지 못한 것이다.
- 고친 문안:
  > The alignment check fixed before computing did not separate aligned from misaligned series: every series other than the reference fell below it, including the repeat at the reference TE, so each was registered on its own FA and the analysis without separately registered series could not be computed under that check (Supplementary Table S10 (d)).
- §2의 결과 문장은 이 문장 바로 뒤에 붙인다.

### 1c. 4절 첫머리: 원뿔 측정의 허용 폭과, 방향과 기울기의 구분

- 3.4절 끝과 4절 첫머리의 «pointed within the algebraic cone»을 «pointed within ±1° of the algebraic cone»으로 바꾼다(정의는 보충 S-M3).
- «…which it did not measure.» 바로 뒤에 넣는다:
  > The measure counts a departure common to all tissue states, which changes intercepts only, as well as one that differs between them, which can change slopes; the measures below concern the second.
- 탐색 문장 괄호에 proportional bias를 넣는다. 단, **먼저 A가 정의를 보고한다.** 무엇을 무엇에 회귀한 기울기인지, 어떤 조직 집합과 단위인지 밝힌다(보충 표 S17(a)에는 정의가 없다).
  - 조직 상태에 따라 이동이 달라지는 정도를 재는 양이 맞으면:
    > (share 1.000; proportional bias −0.023 to +0.035; lowest rank correlation 0.997, no rank reversals)
  - 아니면 넣지 않고, §3의 결과로 대신한다.
- 보충 표 S17(a)에 proportional bias의 정의를 한 줄 덧붙인다.

### 1d. 이론 문장 둘

- **2.2절 «Equivalence» 단락 끝**(«…assessed by simulation (Sections 3.4 and 4).» 뒤)에 넣는다:
  > Models without a tortuosity constraint, such as the Standard Model with free compartmental diffusivities, meet the condition of Proposition 1, so for them the translation is exact wherever the fit is unique, which degeneracies limit (Novikov et al., 2019).
  - 근거: Novikov et al.이 "Standard Model"이라는 이름, 비교환 Gaussian 구획의 합이라는 정의, 퇴화 문제를 다루는 것을 검토 세션이 PMC 전문에서 확인했다(PMC6481929).
  - 참고문헌 추가: Novikov, D.S., Fieremans, E., Jespersen, S.N., Kiselev, V.G., 2019. Quantifying brain microstructure with diffusion MRI: Theory and parameter estimation. NMR in Biomedicine 32, e3998. doi:10.1002/nbm.3998
  - A가 서지 도구로 확인한다. PubMed의 온라인 날짜는 2018-10-15이고, 권호는 32(4)다.
- **6절 결론, 둘째 문장 뒤**에 넣는다:
  > At a single acquisition, a change that scales both tissue weights by one factor, such as a lengthening of T1 shared by both tissue compartments, is exactly equivalent to a change in free water (Proposition 2).
  - 명제 2의 λ 진술을 그대로 옮긴 것이다. 증명된 결과인데 지금 결론에는 없다.
- [저자 결정 필요] 같은 문장을 초록에도 넣을지는 저자가 정한다. 넣으려면 초록에서 같은 낱말 수만큼 줄여야 한다.

## 2. 분석 ①: 다중 TE 자료를 작동하는 정합 검사로 다시 분석 [검토자 제안]

**이유:** 4.3절의 숫자 전부와 5.1절·Box 1의 일부가 각 계열을 따로 정합한 경로에서 나왔다. 해당 숫자는 선형 비율 0.815/0.901, 25.4°, ODI 6–9배, τ·β·ΔTE* 17 ms, E4–E6이다. 원래 고정한 경로는 제공자가 정렬한 계열에 기준 FA의 정합 하나를 똑같이 쓰는 것이었다. 그 경로가 막힌 이유는 자료가 아니라 검사였다.

### 2.1 원인 확인 (서술, 판독 없음)

- 고정 검사의 지표(b0로 정규화한 평균 신호의 복셀 상관)를 다시 계산한다.
- b0가 작은 복셀(예: 마스크 안 b0 하위 1 %)을 뺐을 때 상관이 어떻게 바뀌는지 보고한다.
- 변경 기록에 들어갈 사실 확인만 한다. 새 기준을 정하는 데는 쓰지 않는다.

### 2.2 새 정합 검사 (계산 전에 고정)

- **대상:** 제공자 정렬본의 각 계열(참가자마다 기준 계열을 뺀 9개), 기준 계열(첫 62 ms)의 native 공간.
- **지표:** 각 계열의 평균 b0를 기준 계열의 평균 b0에 FLIRT 6자유도(강체)로 맞춘다. 비용 함수는 normalized mutual information이다. 추정 변환이 뇌 마스크 안에서 만드는 최대 변위(mm)를 잰다.
- **기준:** 최대 변위가 계열 복셀 크기의 절반 이하이면 "정렬됨"이다.
- **실패 대조 (먼저 실행):**
  - (i) 기준 계열과 그 자신: 최대 변위가 0.1복셀 이하여야 한다.
  - (ii) 기준 평균 b0를 x축으로 정확히 2복셀 옮긴 사본: 추정 변위가 2 ± 0.25복셀이어야 하고, "정렬 안 됨"으로 분류돼야 한다.
  - 둘 중 하나라도 어긋나면 검사가 무효다. 멈추고 보고한다.
- **보고만 하는 것:** 같은 TE(62 ms)의 반복 계열 결과. 기대는 "정렬됨"이지만 판독 조건은 아니다.

### 2.3 계산할 경로 (둘 다 계산하고, 결과와 상관없이 둘 다 보고)

- **경로 B (원래 고정한 경로):** 제공자 정렬본 가운데 새 검사를 통과한 계열만 쓴다. 기준 계열의 FA→표준 공간 정합 하나를 모든 계열에 적용한다(applywarp, trilinear). 유효 마스크 규칙(가중치 0.99)과 ROI 규칙은 S-M6 그대로다. 통과한 TE가 선형 판독에 부족하면(고정 규칙의 최소 TE 수) "계산할 수 없음"으로 보고한다.
- **경로 C:** 모든 계열을 2.2의 강체 변환으로 기준 계열에 다시 맞춘 뒤, 기준 계열의 정합 하나를 적용한다.
- **기존 경로 A** (각 계열을 따로 정합): 지금 원고의 값이다. 비교용으로 함께 싣는다.
- 반복 계열(62, 132 ms)도 같은 경로로 처리하고, 재검사 SD는 경로마다 다시 계산한다.
- 경로 사이에 달라지는 것은 정합 단계뿐이다. 코드 차이(diff)를 보고한다.

### 2.4 다시 계산할 양 (경로 A·B·C 나란히)

- 선형 비율 x, y: 고정 규칙(S-M7.4), 전체와 참가자별
- ODI의 TE 범위 ÷ 반복 차이의 중앙값: 참가자별(보충 표 S10)
- TE 기울기, 문헌 예측 대비 기울기 비, 단위 궤적 방향의 중앙값(S-M8.3)
- τ, β, ΔTE*와 그 구간(S-M10 E1), E4, E5, E6(S-M10.5)
- 그림 5와 같은 그림: 경로 B·C(보충 그림)

### 2.5 결과별 문장 (계산 전에 고정)

- **ODI:** 참가자마다 "경로 B의 ODI 범위 비 ÷ 경로 A의 비"를 구한다.
  - 세 참가자 모두 0.5 이하이면:
    > With one registration for all series, the range of ODI across TE was X–Y times the median retest difference, against 6–9 with separately registered series (post hoc), so registration accounted for most of it.
  - 세 참가자 모두 0.8 이상이면 같은 문장 끝을 «…so separate registration did not account for it.»로 한다.
  - 그 밖이면 «…accounted for part of it.»로 한다.
  - 경로 B를 계산할 수 없으면 경로 C로 같은 규칙을 적용하고, 경로 이름을 밝힌다.
- **선형 비율:** 고정 결과(0.815, 0.901)는 그대로 둔다. 4.3절에 경로 B(또는 C)의 값을 «(post hoc, one registration for all series: …)»로 덧붙인다.
- **ΔTE\*:** 경로 B(또는 C)의 값이 경로 A의 90 % 구간(15.1–20.7 ms) 밖이면, 4.3절과 Box 1에 두 값을 함께 쓴다(«about 17 ms with separate registration and about Z ms with one registration»). 구간 안이면 본문은 그대로 두고 보충에만 싣는다.
- **E4–E6과 25.4°:** 보충에 나란히 싣는다. 본문 숫자는 바꾸지 않는다. 경로 B(또는 C)에서 E4 차이의 부호가 바뀌면 4.3절에 한 구절로 밝힌다.

## 3. 분석 ②: 적합기 출력에서 조직 상태에 따라 이동이 달라지는가 (분석 1) [검토자 제안]

**이유:** 측정 불변성(metric invariance)은 "이동이 조직 상태와 무관하다"는 뜻이다. 지금 본문에는 방향(원뿔), 공통 이동 뒤의 크기, 순위만 있다. 심사자가 가장 먼저 물을 질문에 직접 답한다.

### 3.1 입력

- 저장된 stage-3 적합값(S1a): 기준 조직 210개, 이완 설정 17개, 생성기 셋, 쌍 종류 넷(TE만, TR만, 방향 scheme만, 모두 다름)
- 입력 파일의 sha256을 기록한다.

### 3.2 계산

한 칸은 (이완 설정, 생성기, 쌍)이다. 한 칸 안의 조직은 이완을 공유한다.

- **(a) 이동의 조직 의존:**
  - 계산: dᵢ = x_fit,B,i − x_fit,A,i를 참값 (x, y, ODI)에 최소제곱 회귀한다.
  - 보고: 계수와 R², 그리고 Δ_max(조직 범위 안에서 회귀가 예측하는 이동의 최댓값 − 최솟값).
  - 판정 양: Δ_max ÷ (1.96·√2·s), s = 0.031(S1a와 같은 OXF1PRI 재검사 SD 중앙값).
  - 같은 계산을 y 좌표에도 한다. 이때 s는 S-M10에서 계산한 y의 OXF1PRI 재검사 SD 중앙값이다. 평행 이동 주장은 (x, y) 둘 다에 대한 것이다.
- **(b) raw NDI와 대비:**
  - 계산: d_raw = ν_fit,B − ν_fit,A를 참값 (ν, FISO, ODI)에 회귀한다.
  - 판정 양: Δ_max ÷ (1.96·√2·s_raw), s_raw는 S-M10의 raw NDI 재검사 SD다.
  - 같은 칸에서 신호 모형(적합기 없이 식 3)만으로 계산한 d_raw도 함께 보고한다. 명제 3이 예측하는 몫과 적합기가 더한 몫을 나누기 위해서다.
- **(c) 획득 조건별 기울기:**
  - 계산: 각 획득 조건 a에서 x_fit,a를 참값 x에 회귀해 α_a, β_a와 최대 잔차를 구한다. 쌍마다 β_B/β_A를 보고한다.
  - 기준은 β = 1이 아니라 β_B/β_A = 1이다. AMICO는 한 획득 조건 안에서도 β가 1이 아닐 수 있다.
  - 서술이다. 판독은 없다.

### 3.3 실패 대조 (먼저 실행)

- (i) 적합기 대신 명제 1의 대수적 값(이동이 정확히 공통)을 넣는다. 계수가 모두 |·| < 10⁻¹⁰이어야 한다.
- (ii) (i)의 d에 0.10·(x_true − 평균)을 더한다. x 계수가 0.100 ± 10⁻⁶이어야 하고, Δ_max ÷ 폭 > 1이어야 한다.

### 3.4 판독과 결과별 문장 (계산 전에 고정)

- **판독:** 한 쌍 종류에서 모든 칸의 Δ_max ÷ 폭이 1 이하이면 "조직에 따른 이동이 재검사 수준 안"이다. 쌍 종류마다 기준을 넘지 않은 칸의 비율과 최대 비율을 보고한다.
- **모든 쌍 종류가 기준 안이면** 4절 첫머리에 넣는다:
  > At fitter output, the displacement of logit NDI depended on tissue state by at most R of the retest width in any pair class, and that of raw NDI by up to R′ (post hoc; Supplementary Methods S-M10), so differences between tissue states in logit NDI changed by less than retest level.
- **한 쌍 종류라도 넘으면** «…by up to R of the retest width when [pair class]…»로 쓰고, 표 2의 «Associations and ranks» 행 근거 칸에 같은 사실을 적는다.
- 표 2의 같은 행 근거 칸에 logit 대 raw 수치를 한 구절로 넣는다.
- **한계 문장(고정):** «noise-free simulation; not a test of metric invariance in people».

## 4. 분석 ③: 음성 대조 — 참가자별 비율, 동등성 경계, y 실패의 진단 (분석 2) [검토자 제안]

### 4.1 참가자별·중심별 비율 (서술, 판독 없음)

- **참가자별:** 고정 규칙(t₀.₉₇₅(6)·√2·s)으로 참가자마다 1,215단위(27 ROI × 45 중심 쌍)의 재검사 수준 비율을 낸다. 좌표는 x, y, ODI, 그리고 고정 민감도 집합의 raw NDI와 logit FISO다.
- **중심별:** 각 중심이 들어간 9쌍의 단위(729개)로 같은 비율을 낸다.
- **본문:** 4.1절에 «per participant, x a–b, ODI c–d, y e–f (post hoc; Supplementary Table S14)» 한 구절을 넣는다.

### 4.2 동등성 경계 (해석적, post hoc)

- **모형:** 모든 단위에 같은 중심 효과 δ가 있고, 차이의 잡음은 N(δ, 2σ²), s는 자유도 6으로 추정한다고 둔다.
  - 기대 비율은 P(|T₆(λ)| ≤ t₀.₉₇₅(6))이고, λ = δ/(√2σ)다.
  - 기대 비율이 0.90이 되는 δ와, 관측 비율에 해당하는 δ를 구한다.
- **검토 세션 계산값** (A가 같은 식으로 다시 계산해 일치를 확인한다):

  | 비율 | λ | δ/s |
  |---|---|---|
  | 0.90(기준) | 0.768 | 1.086 |
  | 0.931(x) | 0.476 | 0.673 |
  | 0.934(ODI) | 0.437 | 0.618 |
  | 0.844(y) | 1.110 | 1.570 |

  - 효과가 없을 때의 기대 비율은 0.950이다.
- **환산:** 좌표별 c10 재검사 SD의 중앙값으로 logit 단위로 바꾼다. x는 NDI 0.5와 ROI 중앙 NDI에서 NDI 단위로도 바꾼다(dν = ν(1−ν)dx).
- **실패 대조 (모의, 3,645단위, 반복 1,000회):**
  - δ = 0에서 비율의 평균이 0.95 ± 0.005여야 한다.
  - δ = 1.086 s에서 0.90 ± 0.005여야 한다.
  - δ = 2 s에서 0.80보다 작아야 한다(해석값 0.777).
  - 모의에서 s는 ROI 단위마다 한 번 뽑아(자유도 6) 그 ROI의 135단위가 함께 쓴다. 실제 자료와 같은 구조다.
- **본문(3.4절 «Basis of the thresholds» 끝):**
  > Under a center effect common to all units, the expected share falls to 0.90 at about 1.1 retest SDs, and the observed shares correspond to about 0.7 (x) and 0.6 (ODI) retest SDs, or about Q in NDI (post hoc; Supplementary Methods S-M7.3); effects that differ between units are not covered by this bound.
  - 숫자는 A의 재계산값으로 넣는다.

### 4.3 y 실패의 진단 (서술)

- ROI 단위마다 Tong 자료의 평균 FISO와 y의 재검사 수준 비율을 짝짓는다. 다음을 보고한다.
  - Kendall τb
  - 평균 FISO가 중앙값 아래인 단위들과 위인 단위들의 비율
  - x와 logit FISO의 같은 양(비교용)
- **결과별 문장(고정):** y에서 "아래" 비율이 "위" 비율보다 0.05 이상 낮을 때만 4.1절에 넣는다.
  > The units below retest level in y were concentrated where FISO was lowest, where the logit enlarges small absolute differences (post hoc; Supplementary Table S14).
  - 그렇지 않으면 보충에만 싣는다.

## 5. 분석 ④: 다중 TE 잡음 바닥 시뮬레이션 [검토자 제안]

**이유:** 5.3절은 "TE에 따라 SNR이 떨어지고, 잡음 바닥이 궤적에 기여했을 수 있으나 분리하지 않았다"고 쓴다. 적합기만으로 나오는 τ는 이미 있다(0.00035/ms, 관측 0.0029/ms). 같은 틀에 잡음 바닥을 더한다.

### 5.1 입력과 설정

- **환경:** rules_env
- **조직 상태:** S1b의 ROI 대응 조직 상태(다중 TE, TE 62 ms의 관측 ROI 값을 중심 이완에서 식 3으로 참 분율에 대응시킨 것)
- **획득:** 여덟 TE, TR 5.8 s, 방향 scheme은 S-M10.1과 같은 가정
- **SNR:** 참가자별·TE별 측정 세션 SNR(S-M6 방식)
  - TE별 값이 기록돼 있지 않으면, 이미 처리한 b0 파일에서 S-M6 방식으로 계산한다. 새 자료는 읽지 않는다.
- **잡음:** 비중심 χ, L ∈ {1, 4, 16}. 실효 L을 모르므로 셋 다 계산한다.
- **규모:** 단위마다 복셀 200개를 따로 합성하고 적합해 평균한다(ROI 평균과 같은 방식). Monte Carlo 반복은 50회다.
- **양:** 단위별 x·y의 TE 기울기, τ_sim(E1과 같은 잡음 분산 보정), β_sim, 고정 선형 규칙의 비율, ODI 범위 비, 궤적 방향의 중앙값

### 5.2 실패 대조 (먼저 실행)

- (i) SNR = ∞: τ가 S1b의 적합기 단독 값(0.00035/ms)을 ±10 % 안에서 재현해야 한다.
- (ii) SNR = 5: 단위 평균 NDI 또는 FISO가 재검사 폭보다 크게 움직여야 한다. 그래야 이 경로가 잡음 바닥 효과를 만들어 낼 수 있다는 것이 보인다.

### 5.3 판독 (계산 전에 고정한 관례)

- **기준값:** τ_obs는 §2의 경로 A·B·C 가운데 가장 작은 관측 τ로 한다. 잡음 바닥에 가장 유리한 쪽이다.
- **판독 기준:** L 셋 가운데 가장 큰 τ_sim을 쓴다.
  - τ_sim/τ_obs ≤ 0.25이면 "잡음 바닥과 적합기가 설명하지 못함"이다.
  - 0.25–0.75이면 "일부 설명"이다.
  - 0.75 이상이면 "대부분 설명할 수 있음"이다.
- **같은 규칙의 적용:**
  - ODI 범위 비에는 같은 규칙을 쓴다.
  - β에는 부호와 크기를 보고하고, β_sim/β_obs에 같은 경계를 쓴다.

### 5.4 5.1절·4.3절 문장 표 (계산 전에 고정)

| §2 결과(경로 B 또는 C의 β) | §5 결과(τ) | 5.1절 문장 |
|---|---|---|
| β > 0이고 구간이 0을 포함하지 않음 | 설명하지 못함 | Between regions of one brain, one common TE shift did not hold, and neither registration nor the noise floor and fitter accounted for this in reanalysis and simulation, leaving region-specific relaxation as the remaining candidate among those examined (exploratory; post hoc). |
| β > 0이고 구간이 0을 포함하지 않음 | 일부 또는 대부분 | Between regions of one brain, one common TE shift did not hold (exploratory); in simulation the noise floor and fitter could account for [part/most] of the between-region variation (post hoc). |
| 구간이 0을 포함 | (상관없음) | With one registration for all series, the dependence of the TE shift on NDI was not distinguishable from zero (post hoc), so one common TE shift was not shown to fail between regions. |

- 셋째 줄이 나오면 E4–E6의 본문 문장과 «about 17 ms»도 해당 경로의 값과 나란히 다시 쓴다. 문안은 결과를 본 뒤 검토 세션이 대조한다.
- 어느 줄이든 4.3절 끝의 «these analyses do not separate them»은 결과에 맞게 바꾼다.

## 6. 하지 않는 것 (이번 범위 밖) [검토자 제안]

- **보정한 이완 값으로 FISO 부호를 다시 예측하기:** 하지 않는다. 8쌍뿐이고 이미 본 결과(unblinded)라, 맞든 틀리든 무게를 실을 수 없다.
- **심사자가 요구할 때 할 것:**
  - 복셀에서 logit을 먼저 취한 뒤 평균한 결과와의 비교
  - 섬유와 B0의 각도로 ROI 사이 TE 기울기 차이를 설명하는 분석
  - 두 번째 적합기로 실제 자료 분석
  - 식 (7)을 사이트 TE에 회귀해 TE = 0의 집단 차이를 추정하는 정밀도

## 7. 보고 요청

- 이 문서를 정의로 커밋한 해시와 시각(모든 계산보다 앞서야 한다)
- §1: 바뀐 문장 전후, 1a의 단정 표현 목록, 1c의 proportional bias 정의, 1d의 서지 확인 결과
- §2–§5마다:
  - 실패 대조 결과(통과 또는 멈춤)
  - 입력 sha256
  - 결과표
  - 결과별 문장 표에서 고른 줄과 넣은 문장
  - 장부 행(번호는 A가 매긴다)
- §2: 2.1의 원인 확인, 2.2의 검사 결과(계열별 최대 변위), 경로 A·B·C 표, 코드 차이
- 본문 낱말 수와 쪽 수의 변화, H1 관문 결과
- 새 PDF(통독판과 투고판)와 sha256

지킬 것은 그대로: "fixed before computing"만 쓸 것, 우리 출력에 restore/remove/harmonize/correct를 쓰지 않을 것, 판단기 점추정을 결과로 쓰지 않을 것, [A 제안]/[검토자 제안]/[저자 결정] 표시를 유지할 것. 저장소 글에는 «저자 결정»을 쓴다.
