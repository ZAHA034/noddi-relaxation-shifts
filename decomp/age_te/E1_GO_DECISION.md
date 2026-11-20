# E1(나이 × 장비 TE) 돌릴지 판정 — 계획 6a 절 · 2026-11-18

**판정: 돌리지 않는다.** 참가자 정보 파일은 열지 않았다. E1 은 원고의 «앞으로의 과제» 로 쓴다.

규칙은 계획 `ANALYSIS_PLAN.md` 6a 절에 있고, 문헌 검색 전에 고정했다(저장소 밖 기록: 2026-11-18 18:22:45 +0900, 계획 sha256 8b19e79e…). 규칙의 요지는 이렇다. ON-Harmony 나이 범위(대략 20–60 세) 백질 T2 의 나이별 변화를 숫자로 보고한 Q1 문헌 값 가운데 가장 큰 것(L_max, ms/10 년)을 너그럽게 환산한다. 그 낙관적 검정력이 0.5 이상이면(L_max 가 약 0.64 ms/10 년 이상) 돌린다. 쓸 수 있는 문헌 값이 없으면 돌리지 않는다.

## 근거

이번 판정은 «문헌 값이 문턱보다 작음» 이 아니라 «숫자로 된 문헌 값을 하나도 읽지 못함» 에서 나왔다. 규칙의 마지막 줄(쓸 수 있는 값이 없으면 돌리지 않는다)이 적용되었다. 읽은 것 가운데 숫자 없이 방향을 말한 것은 다음과 같다.

| 문헌 | 읽은 범위 | 대상 · 나이 | 말하는 것 |
|---|---|---|---|
| Hagiwara et al., 2021, Invest Radiol (doi:10.1097/RLI.0000000000000720) | 본문 전부(PMC7864648) · 보충 자료(T2 곡선의 식, «Supplemental Digital Content 2»)는 열지 못함 | 건강한 성인 114 명, 21–86 세, 3 T 합성 MRI | «Stable T1 and T2 values from the 20s until 60s in both GM and WM and increases from 60s» |
| Papadaki et al., 2019, J Magn Reson Imaging (doi:10.1002/jmri.26707) | 초록만 | 건강한 성인 90 명, 22–81 세, 1.5 T 다중 에코 스핀 에코 | 긴 T2 의 유의한 증가가 «first noted at 60–69 years» |
| Canales-Rodríguez et al., 2021, Neurobiol Aging (doi:10.1016/j.neurobiolaging.2021.06.002) | 초록만 | 건강한 성인 145 명, 18–60 세, 다중 성분 T2 | 세포 안팎 물의 T2 가 나이에 따라 줄며 남성에서 더 가파름 — 크기는 초록에 없음 |
| Knight et al., 2016, Phys Med Biol (doi:10.1088/0031-9155/61/15/5587) | 본문 전부(PMC5390949) | 대조 37 명, 49–87 세 | 백질 T2 가 나이에 따라 늚. 나이 범위가 규칙 밖이고, 기울기 수치는 그림에만 있음 |

두 연구(Hagiwara, Papadaki)는 20–60 세에서 백질 T2 가 거의 변하지 않는다고 말한다. 이는 «돌리지 않는다» 와 같은 쪽이다. Canales-Rodríguez 는 이 범위에서 세포 안팎 물의 T2 가 준다고 말하지만 크기를 읽지 못했다. 이것이 이 판정에 남은 불확실성이다. 그 수치(또는 Hagiwara 의 보충 식, Hasan et al., 2010 의 곡선 계수)를 얻으면 같은 규칙에 넣어 다시 판정할 수 있다. 규칙은 바꾸지 않는다.

## 검색 기록

PubMed(MCP) 검색 5 회, 2026-11-18. 모두 끝까지 돌았다(Stopped: no).

1. `(T2 relaxation OR T2 relaxometry OR T2 mapping) AND white matter AND (age OR aging OR lifespan) AND (healthy adults OR adult lifespan)` — 158 건 중 25 건의 목록을 보고, 15 건의 서지 · 초록을 열었다.
2. `(white matter[Title/Abstract]) AND (T2 relaxation time[Title/Abstract] OR T2 relaxometry[Title/Abstract] OR quantitative T2[Title/Abstract] OR T2 values[Title/Abstract]) AND (age-related[Title/Abstract] OR aging[Title/Abstract] OR lifespan[Title/Abstract]) AND healthy[Title/Abstract]` — 13 건, 모두 서지 · 초록을 열었다.
3. `Bouhrara[Author] AND (relaxometry OR T2) AND (age OR aging OR lifespan)` — 39 건. 목록만 받고 열지 않았다.
4. `(intra- and extracellular water OR long T2 OR "T2 relaxation") AND white matter AND (age OR aging) AND (myelin water fraction) AND cross-sectional AND adults` — 3 건(이미 연 것).
5. `(synthetic MRI OR QRAPMASTER OR MR fingerprinting OR quantitative MRI) AND (T2) AND white matter AND (age-related OR aging) AND (healthy adults OR normal) AND (relaxation times)` — 18 건, 11 건의 초록을 열었다.

본문을 연 것은 PMC5390949(Knight 2016)와 PMC7864648(Hagiwara 2021) 둘이다. PMC8333472(Bartzokis 1994)는 본문이 비어 있었다.

**막힌 곳(웹 가져오기, 같은 날).**
- Hagiwara 2021 출판사 쪽: HTTP 402. 보충 자료를 열지 못했다.
- Canales-Rodríguez 2021: 출판사 연결 쪽이 «Redirecting» 만 돌려주었다.
- Hasan et al., 2010, Magn Reson Med (doi:10.1002/mrm.22515; 15–59 세 130 명, 백질 T2 의 U 자 곡선): HTTP 403.
- Papadaki 2019 · Canales-Rodríguez 2021 · Hasan 2010 은 PMC 사본이 없다.

**찾지 않은 곳.** Scopus, Web of Science, Google Scholar, 학회 초록, 아카이브(규칙상 쓰지 않음).

## 원고에 쓸 것(보충, 본문 한 문장)

연속 공변량 따름정리의 항 d(R2,en − R2,in)/dz 를 실자료로 재는 일은 앞으로의 과제로 적는다. ON-Harmony 와 같은 여행 피험자 설계라면 나이 SD 10 년, 관측 110 에서도 80 % 최소 검출 크기는 밖 T2 로 10 년당 약 2–4 ms 이다(계획 7 절의 어림). 문헌 둘은 20–60 세의 백질 T2 가 거의 변하지 않는다고 말한다. 한 연구는 이 범위에서 세포 안팎 물의 T2 가 준다고 말하지만 크기를 읽지 못했다. 그래서 더 넓은 나이 범위와, 구획 T2 를 함께 재는 설계가 필요하다. 문헌 인용은 원고 참고문헌 규칙(Q1 · 아카이브 금지)을 따른다. 학술지 등급(SJR 공개 내보내기 사본, 2025 년 값)은 Invest Radiol · J Magn Reson Imaging · Neurobiol Aging 이 Q1, Phys Med Biol 이 Q2(2024 년은 Q1)이다. 원고에는 본문을 직접 읽은 Hagiwara et al., 2021 만 인용한다.

## 덧붙임(판정 뒤, 같은 날)

계획 7 절의 나이 요약(31.5 ± 10.5 세, 여 7)은 앞선 작업 기록에서 옮긴 값이다. 이번에 자료 논문 원문으로 다시 확인하지 않았다(검토 세션은 성별을 «10 명 대 10 명» 으로 적었다). 그래서 원고에는 이 요약을 쓰지 않고 «나이 SD 약 10 년을 가정» 으로만 적는다. 판정은 이 값과 무관하다(규칙의 문헌 값을 하나도 읽지 못함).
