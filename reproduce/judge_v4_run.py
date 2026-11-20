"""Run decomp/remeasure/judge_v4.py with the acquisition schemes in data/protocols.

decomp/stage1/config.py is a locked file (its SHA-256 is listed in docs/fixed_products_and_commits.md), and it
fixes the location of the scheme files relative to the parent of the repository. This script changes
that one location at run time, before the unchanged judge code starts. Nothing else differs.

Usage (the commands of decomp/remeasure/judge_v4.py; reproduce/step2_recompute.py --judge calls them):
  python reproduce/judge_v4_run.py snr-table                 session SNR table (step 1 of V4)
  python reproduce/judge_v4_run.py attr-table <snr sha256>   attribute table of the session pairs (step 2)
  python reproduce/judge_v4_run.py test-cover <arm>          coverage of the simulated test data, arm c3 or c7 (step 3)
  python reproduce/judge_v4_run.py run <arm> <attr sha256>   V4 cells of one arm (step 4)
The outputs go to the folders set by REMEASURE_WORK and DECOMP_WORK; each step stops if the SHA-256 of an
earlier output does not match.
"""
import pathlib
import runpy
import sys

if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
    print(__doc__)
    raise SystemExit(0 if len(sys.argv) > 1 else 2)

REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / "decomp" / d) for d in ("stage1", "s4", "remeasure")]
import config as C  # noqa: E402

C.PROTO_DIR = REPO / "data" / "protocols"
sys.argv = [str(REPO / "decomp" / "remeasure" / "judge_v4.py"), *sys.argv[1:]]
runpy.run_path(sys.argv[0], run_name="__main__")
