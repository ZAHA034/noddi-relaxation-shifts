"""age_te — rules for reading the participant-information tsv and the personal-information guard of the output (plan Sections 3 · 6); reads no data, only parses text it is given."""
import csv
import re

AGE_BIN = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*[-–]\s*(\d+(?:\.\d+)?)\s*$")
SEX = {"M": 0.0, "MALE": 0.0, "F": 1.0, "FEMALE": 1.0}


def parse_participants(text):
    """Text of the participant-information tsv → ({number: age (in units of 10 years)}, {number: female = 1}, format summary). Rules of plan Section 3: age — a number as it is, an interval «a-b» its middle, anything else (n/a and the like) drops that person from E1. Sex — reads only M · F · male · female (case ignored) (run.py does not use it since E2 was dropped · plan change record 2)."""
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    cols = list(rows[0].keys())
    age_col = next((c for c in cols if c.strip().lower().startswith("age")), None)
    sex_col = next((c for c in cols if c.strip().lower() in ("sex", "gender")), None)
    age, female = {}, {}
    fmt = dict(age_column=age_col is not None, sex_column=sex_col is not None, age_bins=0, age_missing=0, sex_missing=0)
    for r in rows:
        pid = r["participant_id"].strip()
        if age_col is not None:
            a = (r[age_col] or "").strip()
            m = AGE_BIN.match(a)
            if m:
                age[pid] = (float(m.group(1)) + float(m.group(2))) / 2 / 10.0
                fmt["age_bins"] += 1
            else:
                try:
                    age[pid] = float(a) / 10.0
                except ValueError:
                    fmt["age_missing"] += 1
        if sex_col is not None:
            s = (r[sex_col] or "").strip().upper()
            if s in SEX:
                female[pid] = SEX[s]
            else:
                fmt["sex_missing"] += 1
    return age, female, fmt



def assert_no_ids(text, ids=()):
    """Keeps participant numbers (sub-…) out of the output — stops if one enters (plan Section 6)."""
    if re.search(r"sub-[0-9A-Za-z]", text) or any(i in text for i in ids):
        raise SystemExit("출력에 참가자 번호가 들어감 — 쓰지 않는다")
