"""Implementation checks of d3 v3: noise-free group estimates within 0.05 of true r, and bit-identical results across two separate processes."""
import json, sys, time, pathlib
import numpy as np
from multiprocessing import get_context
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import d3_v3 as V  # noqa: E402
SIG = np.array([0.0128, 0.0135])


def job(a):
    tag, r_true, scale, noisy = a
    D = V.make_D(0, r_true, scale, noisy=noisy)
    t0 = time.time()
    return tag, r_true, V.prof_cont([D], SIG, V.VARIANTS["dpar1.7|t2iso1000"]), round(time.time() - t0, 1)


if __name__ == "__main__":
    print("V-S", V.build_lib())
    with get_context("fork").Pool(4) as p:
        out = p.map(job, [("N", 1.6, 1.0, False), ("N", 2.64, 1.0, False), ("D1", 2.64, 2.0, True), ("D2", 2.64, 2.0, True)])
    vn = [(t, r, s) for tag, t, r, s in out if tag == "N"]
    ok_n = all(abs(r["rhat"] - t) / t <= 0.05 for t, r, _ in vn)
    for t, r, s in vn:
        print(f"V-N 참값 {t}: r̂ {r['rhat']:.4f} · 상대 차 {(r['rhat'] - t) / t:+.4f} · 구간 [{r['lo']:.3f}, {r['hi']:.3f}] · 평가 {r['n_eval']} · 초 {s}")
    d1, d2 = [r for tag, _, r, _ in out if tag == "D1"][0], [r for tag, _, r, _ in out if tag == "D2"][0]
    ok_d = d1 == d2
    print(f"V-D: r̂ {d1['rhat']:.6f} · 구간 [{d1['lo']:.4f}, {d1['hi']:.4f}] · 위 열림 {d1['open_hi']} · 평가 {d1['n_eval']} · 두 실행 같음 {ok_d} · 초 {[s for tag, _, _, s in out if tag.startswith('D')]}")
    print("V-N", "통과" if ok_n else "실패", "· V-D", "통과" if ok_d else "실패")
    (V.OUT / "d3_v3_check.json").write_text(json.dumps(dict(rule_VN="|r̂ − r| / r ≤ 0.05 (참값 1.6 · 2.64 · 잡음 0 · 주 변형)", rule_VD="같은 입력 두 프로세스 · 결과 비트 단위 같음",
                                                           VN=[dict(r_true=t, **r, sec=s) for t, r, s in vn], VD=dict(run1=d1, run2=d2, equal=bool(ok_d)),
                                                           passed=bool(ok_n and ok_d), script_sha256_d3_v3=__import__("hashlib").sha256(pathlib.Path(V.__file__).read_bytes()).hexdigest()),
                                                      ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
