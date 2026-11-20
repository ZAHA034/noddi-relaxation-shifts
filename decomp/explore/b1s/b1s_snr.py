"""Step 1 of b1s: computes SNR per TE from the b0 repeats of the Wong sub-01 8 TE series, checking each file's sha256."""
import datetime, gzip, hashlib, io, json, pathlib, sys
import numpy as np
import nibabel as nib
REPO = pathlib.Path("<REPO>")
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import analysis as A  # noqa: E402
import rconfig as R  # noqa: E402
W = pathlib.Path("<WORK>/external/wongkornchaovalit2025")
V = pathlib.Path("<WORK>/remeasure/wong/sub-01")
OUT = pathlib.Path("<WORK>/remeasure/explore/b1s")
X = pathlib.Path("<WORK>/stage1_work/drafts/x14")
REG = {p: (s, int(b)) for s, b, p in (ln.split(None, 2) for ln in (REPO / "prereg/artifacts/decomp/wongkornchaovalit2025_sha256.txt").read_text(encoding="utf-8").splitlines())}
now = datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def ext_read(rel, what):
    b = (W / rel).read_bytes(); s = hashlib.sha256(b).hexdigest()
    assert (s, len(b)) == REG[rel], f"⛔ 등록 sha256 · 크기와 다름: {rel}"
    with (X / "EXTERNAL_READ_LOG.jsonl").open("a", encoding="utf-8") as f:
        f.write(json.dumps(dict(시각=now, 자료="Wongkornchaovalit 2025", 읽음=str(W / rel), 읽은_것=what, sha256=s, 바이트=len(b), 등록값과_같음=True,
                                무엇="B1-S — TE 별 b0 SNR · Wong 방향(저자 승인 2026-11-02)", 스크립트=__file__), ensure_ascii=False) + "\n")
    return b


def val_read(p, what):
    b = p.read_bytes()
    with (X / "REREAD_LOG.jsonl").open("a", encoding="utf-8") as f:
        f.write(json.dumps(dict(시각=now, 다시읽음=str(p), sha256=hashlib.sha256(b).hexdigest(), 바이트=len(b), 무엇=what, 스크립트=__file__), ensure_ascii=False) + "\n")
    return nib.Nifti1Image.from_bytes(gzip.decompress(b)).get_fdata()


rows = {}
for te in (62, 72, 82, 92, 102, 112, 122, 132):
    pre = f"MTE-dMRI/sub-01/dwi/sub-01_acq-TE{te}_dwi"
    bv = np.loadtxt(io.BytesIO(ext_read(pre + ".bval", "bval 전체"))).ravel()
    img = nib.Nifti1Image.from_bytes(gzip.decompress(ext_read(pre + ".nii.gz", "4D 전체(b0 부피만 씀)")))
    data = np.asarray(img.dataobj)
    b0 = data[..., bv <= R.AMICO_B0_THR].astype(float)
    sd = V if te == 62 else V / f"own_TE{te}"
    msk = val_read(sd / "nodif_brain_mask.nii.gz", "B1-S SNR — 계열 마스크"); fa = val_read(sd / "dti_FA.nii.gz", "B1-S SNR — 계열 FA")
    assert msk.shape == b0.shape[:3] == fa.shape, (te, msk.shape, b0.shape, fa.shape)
    s = A.snr(b0, msk, fa)
    rows[str(te)] = dict(s, b_values=sorted({float(x) for x in np.round(bv)}), n_vol=int(len(bv)))
    if te == 92:
        bvec = np.loadtxt(io.BytesIO(ext_read(pre + ".bvec", "bvec 전체"))).reshape(3, -1).T
        np.savez(OUT / "wong_sub01_TE92_scheme.npz", bvals=bv, bvecs=bvec)
    print(te, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in rows[str(te)].items() if k != "b_values"}, rows[str(te)]["b_values"])
(OUT / "b1s_snr.json").write_text(json.dumps(dict(made=now, script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
                                                  definition="analysis.snr(b0 부피들, 계열 마스크, 계열 dti_FA) — b0 = bval ≤ 10(AMICO_B0_THR)", te=rows), ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
