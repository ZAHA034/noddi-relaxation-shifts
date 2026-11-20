"""b2k step 1 for the b2000 noise floor: reads headers, bval and sidecars of 109 V4 first sessions into a session table, without voxel values."""
import collections, datetime, hashlib, json, os, pathlib, sys

import numpy as np

REPO = pathlib.Path("<REPO>")
sys.path.insert(0, str(REPO / "decomp" / "remeasure"))
import rconfig as R  # noqa: E402
import guard as G  # noqa: E402
import nibabel as nib  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import b2k_core as C  # noqa: E402

X14 = pathlib.Path("<WORK>/stage1_work/drafts/x14")
OUT = R.WORK / "explore" / "b2000"
TAG = "봉인 해제 — 사후 탐색 · 판정 불변 #34"
SHA_ATTR = "5184137854abc69710904889503eaf8dbb111ca9736e46af13b49067ea65b90d"
EXPECT = {"NOT1ACH": 20, "NOT2ING-89": 4, "NOT2ING-98": 16, "NOT4GEP": 10, "OXF1PRI": 20, "OXF2PRI": 20, "OXF3TRI": 10, "OXF4GEP": 9}
SIDE_KEYS = ("Manufacturer", "ManufacturersModelName", "SoftwareVersions", "MagneticFieldStrength", "ReceiveCoilName", "ReceiveCoilActiveElements",
             "CoilCombinationMethod", "ImageType", "ScanOptions", "SequenceVariant", "SequenceName", "PulseSequenceDetails", "PartialFourier",
             "PartialFourierDirection", "ParallelReductionFactorInPlane", "ParallelAcquisitionTechnique", "MultibandAccelerationFactor",
             "NonlinearGradientCorrection", "EchoTime", "RepetitionTime", "PhaseEncodingDirection", "AcquisitionMatrixPE", "ReconMatrixPE",
             "BaseResolution", "PercentPhaseFOV", "SliceThickness", "SpacingBetweenSlices", "EffectiveEchoSpacing", "TotalReadoutTime",
             "PixelBandwidth", "PhilipsRescaleSlope", "PhilipsRescaleIntercept", "PhilipsScaleSlope", "UsePhilipsFloatNotDisplayScaling",
             "ImageComments", "SeriesDescription", "ProtocolName", "ConversionSoftware", "ConversionSoftwareVersion")


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def sha(p: pathlib.Path) -> tuple[str, int]:
    h, n = hashlib.sha256(), 0
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
            n += len(blk)
    return h.hexdigest(), n


def unit_of(sub: str, ses: str) -> str:
    dev = R.device_of(ses)
    return R.not2ing_protocol(sub) if dev == "NOT2ING" else dev


def sessions() -> list[tuple[str, str]]:
    ap = R.WORK / "onharmony" / "v4" / "attr_table.json"
    b = ap.read_bytes()
    assert hashlib.sha256(b).hexdigest() == SHA_ATTR, "attr_table.json sha256 가 기록과 다르다"
    with (X14 / "REREAD_LOG.jsonl").open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(dict(시각=now(), 다시읽음=str(ap), sha256=SHA_ATTR, 바이트=len(b), 무엇="b2000 단계 1 — 첫 세션 목록(§27)",
                                 스크립트=__file__), ensure_ascii=False) + "\n")
    at = json.loads(b)
    sess = sorted({(sub, s) for d in at["pairs"].values() for sub, x in d.items() for s in x["sessions"]})
    got = collections.Counter(unit_of(a, b) for a, b in sess)
    assert dict(got) == EXPECT, f"세션 수가 §27 과 다르다: {dict(got)}"
    return sess


def jsonable(v):
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return None if not np.isfinite(v) else float(v)
    if isinstance(v, np.ndarray):
        return [jsonable(x) for x in v.tolist()]
    if isinstance(v, (list, tuple)):
        return [jsonable(x) for x in v]
    return v


def main():
    basis = G.require()
    sess = sessions()
    rows = {}
    for sub, ses in sess:
        f_nii = R.ONH_SEALED / sub / ses / "dwi" / f"{sub}_{ses}_dir-AP_dwi.nii.gz"
        f_bval = R.ONH_META / sub / ses / "dwi" / f"{sub}_{ses}_dir-AP_dwi.bval"
        f_json = R.ONH_META / sub / ses / "dwi" / f"{sub}_{ses}_dir-AP_dwi.json"
        logs = []
        for f in (f_nii, f_bval, f_json):
            h, n = sha(f)
            logs.append(dict(시각=now(), utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"), 세션=f"{sub}_{ses}",
                             파일=os.path.relpath(f, REPO), sha256=h, 바이트=n,
                             무엇=f"{TAG} · b2000 단계 1 헤더(NIfTI 는 헤더만 · 화소 값 읽지 않음)", 근거=basis, 스크립트=__file__))
        with (X14 / "ONH_NEW_READ_LOG.jsonl").open("a", encoding="utf-8") as fh:
            for e in logs:
                fh.write(json.dumps(e, ensure_ascii=False) + "\n")
        bval = np.loadtxt(f_bval).ravel()
        i0, i2 = C.shells(bval)
        img = nib.load(str(f_nii))
        hd = img.header
        slope, inter = hd.get_slope_inter()
        meta = json.loads(f_json.read_bytes())
        assert img.shape[3] == len(bval), (sub, ses, img.shape, len(bval))
        rows[f"{sub}_{ses}"] = dict(
            unit=unit_of(sub, ses), device=R.device_of(ses), n_vol=int(len(bval)), n_b0=int(len(i0)), n_b2k=int(len(i2)),
            i_b0=i0.tolist(), i_b2k=i2.tolist(), b_values=sorted({int(round(b)) for b in bval}),
            b2k_values=sorted({int(round(b)) for b in bval[i2]}),
            shape=[int(x) for x in img.shape], zooms=[round(float(z), 4) for z in hd.get_zooms()],
            datatype=str(hd.get_data_dtype()), scl_slope=jsonable(slope), scl_inter=jsonable(inter),
            sform_code=int(hd["sform_code"]), qform_code=int(hd["qform_code"]), affine=[[round(float(v), 5) for v in r] for r in img.affine],
            side={k: jsonable(meta.get(k, "칸 없음")) for k in SIDE_KEYS}, side_keys=sorted(meta),
            sha256={"nii": logs[0]["sha256"], "bval": logs[1]["sha256"], "json": logs[2]["sha256"]})
        del meta, img
    summ = {}
    for u in EXPECT:
        vs = [v for v in rows.values() if v["unit"] == u]
        keys = ("n_vol", "n_b0", "n_b2k", "b2k_values", "shape", "zooms", "datatype", "scl_slope", "scl_inter")
        s = {k: dict(collections.Counter(json.dumps(v[k]) for v in vs)) for k in keys}
        s["side"] = {k: dict(collections.Counter(json.dumps(v["side"][k], ensure_ascii=False) for v in vs)) for k in SIDE_KEYS}
        s["n_sessions"] = len(vs)
        summ[u] = s
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "headers.json").write_text(json.dumps(dict(made=now(), plan="§27 (2026-11-03 21:23:18 KST)", tag=TAG,
                                                      script_sha256=sha(pathlib.Path(__file__))[0], core_sha256=sha(HERE / "b2k_core.py")[0],
                                                      basis=basis, n_sessions=len(rows), by_unit=summ, sessions=rows),
                                                 ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    for u, s in summ.items():
        print(u, s["n_sessions"], "vol", s["n_vol"], "b0", s["n_b0"], "b2k", s["n_b2k"], "shape", s["shape"], "zooms", s["zooms"], s["datatype"])


if __name__ == "__main__":
    main()
