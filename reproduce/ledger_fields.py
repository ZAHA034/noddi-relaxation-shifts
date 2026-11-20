"""Field readers of the ledger: json:, jsonl:, expr:, md:, lit:, none:. Used by step1_check_numbers.py."""
from __future__ import annotations

import csv, hashlib, json, math, pathlib, re, statistics, sys
from decimal import Decimal, ROUND_HALF_UP

HERE = pathlib.Path(__file__).resolve().parent
TAG = re.compile(r"\[(N-\d{3,4})\]")
FRAC = "¼½¾⅐⅑⅒⅓⅔⅕⅖⅗⅘⅙⅚⅛⅜⅝⅞"
NUM = re.compile(r"(?<![\w@./-])[−+±-]?(?:\d+(?:[.,]\d+)*(?:\s?%|°)?(?![\w])|[" + FRAC + r"])")   # Fraction characters (½ ⅓ …) are also numbers
NUMTOK = re.compile(r"\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")
SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789-+")
MS_WORDS = {"no": 0, "none": 0, "single": 1, "once": 1, "twice": 2}      # Only on the manuscript-notation side
STRICT = re.compile(r"<!--\s*ledger:strict\s*-->(.*?)<!--\s*/ledger:strict\s*-->", re.S)
WORDS = {"zero": 0, "one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10,
         "eleven": 11, "twelve": 12, "twenty": 20, "half": 0.5,
         # Korean number words: they read numbers written in Korean in the cited working records
         "하나": 1, "한": 1, "둘": 2, "두": 2, "셋": 3, "세": 3, "넷": 4, "네": 4,
         "다섯": 5, "여섯": 6, "일곱": 7, "여덟": 8, "아홉": 9, "열": 10, "절반": 0.5}
EXCEPTIONS = [
    re.compile(r"<!--.*?-->", re.S),                          # Annotation (including the strict-zone marker)
    re.compile(r"\$\$.*?\$\$", re.S),                          # Formula (display) — constants inside a formula are not ledger numbers
    re.compile(r"(?<![\$\\])\$[^$\n]+\$(?!\$)"),                # Formula (in text)
    re.compile(r"`[^`\n]*`"),                                 # Code · file path
    re.compile(r"https?://\S+|doi:\S+|\b10\.\d{4,}/\S+"),     # URL · DOI
    re.compile(r"\[src:[^\]]*\]\{\.src\}"),                   # Within the basis marker
    re.compile(r"\[@[^\]]*\]|-?@[A-Za-z][\w:.-]*"),            # Citation
    re.compile(r"\{#[^}]*\}"),                                # Cross-reference label
    re.compile(r"\b(?:Sections?|Fig(?:ure)?s?\.?|Tables?|Box(?:es)?|Appendix|Propositions?|§)\s*S?[\d.]+[a-z]?(?:\s*(?:[–-]|and|,)\s*S?[\d.]+[a-z]?)*"),
    re.compile(r"\b(?:Eqs?\.|Equations?)\s*\(?\d+\)?(?:\s*(?:[–-]|and|,)\s*\(?\d+\)?)*"),
    re.compile(r"^#+\s.*$", re.M),                             # Title line
    re.compile(r"^(?:\s*>)*\s*\d+\.(?=\s)", re.M),                # Leading number of a numbered list (1. 2. …)
    re.compile(r"⟦[^⟧]*⟧"),
    re.compile(r"\[(?:검토자 제안|A 고침|A 제안|저자 결정)[^\]]*\]"),
    re.compile(r"\bb\s?=\s?(?:0|1000|2000|3000)(?:\s?(?:and|or|,)\s?(?:1000|2000|3000))*(?:\s?s/mm²)?"),   # Shell name  — before the year form (2000)
    re.compile(r"\b\d{4}-\d{2}-\d{2}\b|\b\d{2}:\d{2}(?::\d{2})?\b"),   # Date · time
    re.compile(r"\b(?:19|20)\d{2}[a-z]?\b"),                   # Year
    re.compile(r"\b[0-9a-f]{7,40}\b"),                         # Commit hash
    re.compile(r"\b\d{2}-to-\d{2}-channel\b|\b\d{2}-channel\b"),   # Coil names (32-channel · 64-channel · 64-to-32-channel)
    re.compile(r"√\d+"),                                       # Formula constant (√2)
    re.compile(r"\b(?:1\.5|3|7)\s?T\b"),                       # Field strength
    re.compile(r"\b9[05]\s?%(?=\s+(?:[\w-]+\s+){0,2}(?:intervals?|CIs?)\b)"),   # Interval level
    re.compile(r"\bClaude Opus \d+(?:\.\d+)?\b|\bclaude-opus-\d+(?:-\d+)*\b"),   # Model name · version (tool and version notation of the v12 · Elsevier GenAI policy — not a result number)
    TAG,
]
NUMWORD = re.compile(r"\b(?:two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|thirteen|fourteen|fifteen|sixteen|seventeen|eighteen|"
                     r"nineteen|twenty|thirty|forty|fifty|hundred|thousand)\b", re.I)   # Count words of the strict zone (author decision E15 · 2026-11-04)
WORD_EXCEPTIONS = [   # Number of sentence structures — not the number of data (scanners · participants · cells · sessions) but the composition of text · model · design (E15 · only countable data numbers are ledger rows)
    (re.compile(r"\bthe two\b", re.I), "refers to two items named just before"),
    (re.compile(r"\b(?:three|four)\s+(?:contributions|questions|items|checks|results stand)\b", re.I), "structure of the text (contributions, questions, list items, checks; «Four results stand» in 5.1 — R35, S84 5 (2))"),
    (re.compile(r"\btwo (?:features|kinds)\b", re.I), "structure of the text (Introduction of v11: «Two features …», «Two kinds of tools …»)"),
    (re.compile(r"\bthree fractions\b|\btwo log-ratios\b|\ball three compartmental\b", re.I), "structure of the model (three fractions, two log-ratios, three compartments)"),
    (re.compile(r"\btwo (?:groups|acquisitions)\b", re.I), "general wording (any two groups, two acquisitions)"),
    (re.compile(r"\b(?:two|three) (?:mechanisms|noise assumptions|further estimators)\b", re.I), "structure of the discussion"),
    (re.compile(r"\btwo hemispheric indicators\b", re.I), "design of a cited study"),
    (re.compile(r"\bwith two candidate combinations\b|\bcells with four\b", re.I), "number of candidate coil combinations (2, 4 = 2¹, 2²)"),
]
FORBIDDEN = [
    (re.compile(r"\b(?:restor|remov|harmoni[sz]|correct)\w*", re.I), "forms of restore, remove, harmonize, correct (not used for this study's outputs; kept when naming another method)"),
    (re.compile(r"\bsimulation-based (?:identifiability )?judge", re.I), "write simulation-trained judge"),
    (re.compile(r"\bcertif\w*", re.I), "forms of certify"),
    (re.compile(r"\b(?:pre-?registered|pre-specified|a priori)\b", re.I), "write fixed before computing (registered applies only to the earlier study)"),
    (re.compile(r"\b(?:travelling|centres?|centred|judgements?|fibres?|grey|behaviour|analys(?:e|ed|ing))\b", re.I), "American spelling"),
    (re.compile(r"\bpost-hoc\b", re.I), "post hoc (no hyphen)"),
    (re.compile(r"\b(?:novel|novelty|for the first time|unprecedented|state-of-the-art)\b", re.I), "novelty words"),
    (re.compile(r"\bFWF\b"), "write FISO"),
]


def norm_text(s: str) -> str:
    """Text for comparison — drops spaces · ** · ` and unifies the forms of the minus sign."""
    return re.sub(r"[\s   ]+", "", s.replace("**", "").replace("`", "")).replace("−", "-")


def numbers_in(s: str, ms: bool = False) -> list[float]:
    """Numbers in the text (thousands commas · signs · the form 10⁻⁴ · English and Korean numerals). Numbers attached to words (OXF1PRI · TE62 · c09) are left out."""
    t = s.replace("−", "-").replace("**", "")
    t = re.sub(r"(?:(\d+(?:\.\d+)?)\s*[×x]\s*)?10([⁻⁺]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+)",
               lambda m: f"{m.group(1) or '1'}e{m.group(2).translate(SUP)}", t)
    out = []
    for m in NUMTOK.finditer(t):
        i = m.start()
        if i > 0 and (t[i - 1].isalpha() or t[i - 1] in "_."):
            continue
        v = float(m.group(0).replace(",", ""))
        if i > 0 and t[i - 1] in "-+" and (i == 1 or not t[i - 2].isalnum()):
            v = -v if t[i - 1] == "-" else v
        out.append(v)
    for ch, v in zip(FRAC, (1/4, 1/2, 3/4, 1/7, 1/9, 1/10, 1/3, 2/3, 1/5, 2/5, 3/5, 4/5, 1/6, 5/6, 1/8, 3/8, 5/8, 7/8)):
        for m in re.finditer(re.escape(ch), t):
            i = m.start()
            out.append(-v if i > 0 and t[i - 1] in "-" else v)
    for w, v in (WORDS | MS_WORDS if ms else WORDS).items():
        if re.search(rf"(?<!\w){re.escape(w)}(?!\w)", t, re.I):
            out.append(float(v))
    return out


def is_number(part: str) -> bool:
    try:
        float(part.replace("−", "-").replace(",", ""))
        return True
    except ValueError:
        return False


def get_path(obj, path: str):
    cur = obj
    for m in re.finditer(r'\["([^"]*)"\]|\[\'([^\']*)\'\]|\[(\d+)\]|([^.\[\]]+)', path):
        q1, q2, idx, key = m.groups()
        if idx is not None:
            cur = cur[int(idx)]
        else:
            cur = cur[q1 if q1 is not None else q2 if q2 is not None else key]
    return cur


def load_src(src_dir: pathlib.Path, rel: str, sha: str = ""):
    base = src_dir / rel.replace("../", "__up__/")
    p = base / sha[:12] if sha and (base / sha[:12]).is_file() else base
    return p, p.read_bytes()


def split_parts(v: str) -> list[str]:
    return [p.strip() for p in str(v).split("|")]


def _expr(doc, rest: str, coerce: bool):
    def sub(m):
        v = get_path(doc, m.group(1))
        if coerce and isinstance(v, str):
            try:
                v = float(v)
            except ValueError:
                pass
        return repr(v)
    expr = re.sub(r"\{json:([^}]+)\}", sub, rest)
    env = {"__builtins__": {}, "round": round, "abs": abs, "math": math, "sum": sum, "len": len, "min": min, "max": max,
           "sorted": sorted, "float": float, "int": int, "str": str, "any": any, "all": all, "statistics": statistics,
           "True": True, "False": False, "None": None}
    return eval(expr, env, {})


def eval_field(src_dir: pathlib.Path, row: dict):
    """Basis field → recomputed value. For md it also compares the pieces and returns the "original value". For lit · none, None."""
    rel, field, sha = row["근거 파일"], row["근거 필드"], row["근거 파일 sha256"]   # Korean column names of the ledger rows: source file, source field, SHA-256 of the source file
    kind, _, rest = field.partition(":")
    if kind in ("lit", "none"):
        return None
    p, b = load_src(src_dir, rel, sha)
    if kind == "json":
        return get_path(json.loads(b), rest)
    if kind == "jsonl":
        kv, _, path = rest.partition(":")
        k, _, v = kv.partition("=")
        for ln in b.decode("utf-8").splitlines():
            r = json.loads(ln)
            if str(r.get(k)) == v:
                return get_path(r, path)
        raise KeyError(f"no line with {k}={v} in the jsonl file")
    if kind == "expr":
        doc = json.loads(b)
        try:
            return _expr(doc, rest, coerce=False)
        except TypeError:
            return _expr(doc, rest, coerce=True)
    if kind == "md":
        line, _, frag = rest.partition(":")
        lines = b.decode("utf-8").splitlines()
        text = lines[int(line) - 1]
        if frag not in text:
            raise ValueError(f"{rel}:{line} does not contain the fragment «{frag[:60]}»")
        nf, nums = norm_text(frag), numbers_in(frag)
        missing = []
        for part in split_parts(row["원 값"]):
            if not part or norm_text(part) in nf:
                continue
            pn = numbers_in(part)
            if pn and all(any(math.isclose(x, y, rel_tol=1e-9, abs_tol=1e-12) for y in nums) for x in pn):
                continue
            missing.append(part)
        if missing:
            raise ValueError(f"{rel}:{line}: the fragment lacks the parts {missing} of the exact value")
        return row["원 값"]
    raise ValueError(f"unknown form of source field: {field}")


def _close(a, b) -> bool:
    if isinstance(a, bool) or isinstance(b, bool) or a is None or b is None:
        return a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_close(a[k], b[k]) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_close(x, y) for x, y in zip(a, b))
    return a == b


def same_value(src, ledger: str) -> bool:
    s = str(ledger).strip()
    if isinstance(src, (dict, list)) and s[:1] in "{[":
        try:
            return _close(src, json.loads(s))
        except ValueError:
            pass
    if isinstance(src, (list, tuple)):
        parts = split_parts(s)
        return len(parts) == len(src) and all(same_value(x, y) for x, y in zip(src, parts))
    if isinstance(src, bool) or src is None:
        return {"true": True, "false": False, "null": None, "none": None}.get(s.lower(), "?") == src
    try:
        return math.isclose(float(src), float(s.replace("−", "-").replace(",", "")), rel_tol=1e-9, abs_tol=1e-12)
    except (TypeError, ValueError):
        return str(src).strip() == s


