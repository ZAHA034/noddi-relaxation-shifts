"""Tong 판정 — 고정 코드 decomp/remeasure/judge_tong.py 의 명령줄(:223–229)이 의도한 계산을 밖에서 그대로 부른다(저자가 (가)를 고름 · 2026-10-31 · 이탈로 기록).
judge_tong.py 는 고치지 않는다(가드 목록 32 파일). 명령줄과 다른 점은 하나 — :227 의 load 반환을 셋으로 받는다
(load 는 08a189c 부터 (values, rep, missing) 을 돌려주는데 :227 은 둘로 받아 ValueError · 2026-10-30 18:09 종료 1 · 결과 0).
결과 JSON 은 값 영역 WORK/tong/judge/tong_judge.json 에 한 번만(있으면 멈춤) · sha256 을 열람 기록(external «Tong 2020» «열람(새 분석)»)에
run.log 로 적는다(판정 코드와 같은 열람 기록 잠금 · 저장소 쓰기 잠금은 도는 V4 실행기가 쥐고 있다). 표준 출력에는 경로 · sha256 만.
    cd dmri_cov && CUDA_VISIBLE_DEVICES= uv run --frozen python <이 파일>"""
import hashlib
import json
import pathlib
import sys

REPO = pathlib.Path("<repo>")
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import guard as G  # noqa: E402
import judge_tong as JT  # noqa: E402
import rconfig as R  # noqa: E402
import run as RUN  # noqa: E402

work = R.WORK
dst = work / "tong" / "judge" / "tong_judge.json"
if dst.exists():
    raise SystemExit(f"⛔ {dst} 가 이미 있다 — 다시 쓰지 않는다")
# ── judge_tong.py :224 – 229 그대로(:227 만 셋으로 받음) ─────────────────────────
v, rp, miss = JT.load(work)
res = dict(main=JT.judge(v, rp), without_c09=JT.judge(v, rp, tuple(c for c in JT.CENTERS if c != "c09")), missing_sessions=miss)
vs, rps, _miss_all = JT.load(work, fname="roi_allshell.json")
res["all_shells_descriptive"] = JT.judge(vs, rps)
text = json.dumps(res, ensure_ascii=False, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
# ── 한 번만 저장 · 열람 기록 ─────────────────────────────────────────────────────
dst.parent.mkdir(parents=True, exist_ok=True)
with open(dst, "x", encoding="utf-8") as fh:
    fh.write(text + "\n")
RUN.log("tong", "판정 — judge_tong 결과(tong_judge.json) · 고정 코드 judge_tong.py:227 명령줄 결함 때문에 같은 함수(load · judge)를 밖에서 부름 — "
        "이탈(저자가 (가)를 고름 · 2026-10-31) · 계산은 명령줄과 같음", [dst], G.require(), "판정")
print(dst, hashlib.sha256(dst.read_bytes()).hexdigest())
