"""Trains the s4 models on CPU: C3 and C7 ensembles, oracle NPE and observation density q(x | a), with sbi 0.27.0 and torch."""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s4cfg as S  # noqa: E402

torch.set_num_threads(int(os.environ.get("S4_THREADS", "16")))

ARCH = dict(hidden_features=128, num_transforms=5, emb_layers=3, emb_hidden=128, lr=5e-4, batch=256, val=0.1, patience=20, max_epochs=500)
_SPEED = S.RESULTS / "speed.json"
if _SPEED.exists():                                    # SPEC_S4 §10-1 — epoch upper limit fixed after the speed measurement (before the results)
    import json as _json
    ARCH["max_epochs"] = min(500, int(_json.loads(_SPEED.read_text(encoding="utf-8"))["estimate"]["epoch_cap"]))
SEEDS = (11, 12, 13, 14, 15)
MODELS = S.WORK / "models"


def prior_box(lo, hi):
    from sbi.utils import BoxUniform
    return BoxUniform(torch.as_tensor(lo, dtype=torch.float32), torch.as_tensor(hi, dtype=torch.float32))


def theta_box():
    lo = np.r_[-S.D_TRAIN, [S.BOX["NDI"][0], S.BOX["FISO"][0], S.BOX["ODI"][0]]]
    hi = np.r_[S.D_TRAIN, [S.BOX["NDI"][1], S.BOX["FISO"][1], S.BOX["ODI"][1]]]
    return lo, hi


def inputs(sim, model):
    cond = sim["cond"].astype(np.float32)
    if model == "c3":
        return np.c_[sim["x_c3"], cond]
    if model == "c7":
        return np.c_[sim["x_c7"], cond]
    if model == "oracle":
        return np.c_[sim["oracle"], cond]
    raise ValueError(model)


def build(model, dim_x):
    from sbi.inference import NLE, NPE
    from sbi.neural_nets import likelihood_nn, posterior_nn
    from sbi.neural_nets.embedding_nets import FCEmbedding
    if model == "q":
        de = likelihood_nn(model="nsf", hidden_features=ARCH["hidden_features"], num_transforms=ARCH["num_transforms"],
                           z_score_theta="independent", z_score_x="independent")
        c = np.zeros((1, len(S.COND_NAMES)))
        return NLE(prior=prior_box(c[0] - 1e6, c[0] + 1e6), density_estimator=de, device="cpu")
    emb = FCEmbedding(input_dim=dim_x, output_dim=ARCH["emb_hidden"], num_layers=ARCH["emb_layers"], num_hiddens=ARCH["emb_hidden"])
    de = posterior_nn(model="nsf", hidden_features=ARCH["hidden_features"], num_transforms=ARCH["num_transforms"], embedding_net=emb,
                      z_score_theta="independent", z_score_x="independent")
    return NPE(prior=prior_box(*theta_box()), density_estimator=de, device="cpu")


def fit(model, seed_i, probe=None):
    if torch.cuda.is_available():
        raise SystemExit("⛔ CUDA 가 보인다 — CPU 만(PLAN §9-1 #49). CUDA_VISIBLE_DEVICES= 로 다시 돌려라")
    src = "oracle_train" if model == "oracle" else "train"
    sim = dict(np.load(S.WORK / f"sim_{src}.npz", allow_pickle=False))
    n = len(sim["theta"])
    rows = np.arange(n) if probe is None else np.random.default_rng(0).choice(n, min(probe[0], n), replace=False)
    torch.manual_seed(SEEDS[seed_i])
    np.random.seed(SEEDS[seed_i])
    if model == "q":
        th = torch.as_tensor(sim["cond"][rows], dtype=torch.float32)
        x = torch.as_tensor(sim["x_c3"][rows], dtype=torch.float32)
        inf = build(model, x.shape[1])
    else:
        x = torch.as_tensor(inputs(sim, model)[rows], dtype=torch.float32)
        th = torch.as_tensor(sim["theta"][rows], dtype=torch.float32)
        inf = build(model, x.shape[1])
    inf.append_simulations(th, x)
    t0 = time.time()
    epochs = probe[1] if probe else ARCH["max_epochs"]
    est = inf.train(training_batch_size=ARCH["batch"], learning_rate=ARCH["lr"], validation_fraction=ARCH["val"],
                    stop_after_epochs=(10 ** 9 if probe else ARCH["patience"]), max_num_epochs=epochs, show_train_summary=False)
    sec = time.time() - t0
    n_ep = len(inf._summary["epochs_trained"]) and int(inf._summary["epochs_trained"][-1])
    rec = dict(UTC=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), model=model, seed=SEEDS[seed_i], rows=len(rows), epochs=n_ep,
               seconds=round(sec, 1), sec_per_epoch=round(sec / max(n_ep, 1), 2), threads=torch.get_num_threads(), load1=os.getloadavg()[0],
               probe=bool(probe))
    if not probe:
        MODELS.mkdir(parents=True, exist_ok=True)
        path = MODELS / f"{model}_{seed_i}.pt"
        torch.save(dict(state=est.state_dict(), arch=ARCH, model=model, seed=SEEDS[seed_i], rows=len(rows)), path)
        rec["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        torch.save(inf.build_posterior() if model != "q" else est, MODELS / f"{model}_{seed_i}.obj.pt")
    (S.WORK / "logs").mkdir(parents=True, exist_ok=True)
    with open(S.WORK / "logs" / "train.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(json.dumps({k: v for k, v in rec.items() if k != "sha256"}), flush=True)
    return rec


if __name__ == "__main__":
    a = sys.argv[1:]
    probe = None
    if "--probe" in a:
        k = a.index("--probe")
        probe = (int(a[k + 1]), int(a[k + 2]))
        a = a[:k]
    fit(a[0], int(a[1]) if len(a) > 1 else 0, probe)
