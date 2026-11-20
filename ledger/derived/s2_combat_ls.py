"""Derived ledger values: summary of the ComBat location-and-scale arm (combat_ls) of analysis S2 of the measurement-function environment (no new computation: minimum and maximum of a ratio and the number of cells outside 0.80–1.25)."""
import hashlib, json, pathlib, subprocess, sys

sys.dont_write_bytecode = True
HERE = pathlib.Path(__file__).resolve().parent
S2 = HERE.parents[1] / "env" / "out" / "s2.json"
r = subprocess.run([sys.executable, "-I", "-B", str(HERE.parent / "reread_log.py"), str(S2), "장부 도출 값(s2_combat_ls.py) — S2 combat_ls 팔 요약 · 새 계산 없음"], capture_output=True, text=True)
h = r.stdout.strip().splitlines()[-1] if r.returncode == 0 and r.stdout.strip() else None
raw = S2.read_bytes()
if h != hashlib.sha256(raw).hexdigest():
    raise SystemExit("STOP: s2.json 기록 실패 또는 바뀜")
by = json.loads(raw)["by"]
SC, GEN = ("fig3", "l2en", "l2both"), ("T-in", "T-app", "T-none")
LOC = ("none", "ctrl_offset", "combat_loc", "combat_train_ctrl_loc")
ls = [by[f"pop|x|two|{s}|{g}|combat_ls"]["ratio"] for s in SC for g in GEN]
loc = [by[f"pop|x|two|{s}|{g}|{m}"]["ratio"] for s in SC for g in GEN for m in LOC]
res = dict(source={"s2.json": h}, cells="pop|x|two × (fig3, l2en, l2both) × (T-in, T-app, T-none)",
           combat_ls=dict(ratio_min=min(ls), ratio_max=max(ls), n=len(ls), n_outside_band=sum(1 for v in ls if v < 0.80 or v > 1.25)),
           location_check=dict(ratio_min=min(loc), ratio_max=max(loc), n=len(loc), n_outside_band=sum(1 for v in loc if v < 0.80 or v > 1.25)))
(HERE / "s2_combat_ls.json").write_text(json.dumps(res, indent=1) + "\n", encoding="utf-8")
print(json.dumps(res))
