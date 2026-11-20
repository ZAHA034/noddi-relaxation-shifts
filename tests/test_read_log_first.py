"""The read log is written before the file is read, and a log that cannot be written stops the read."""
from __future__ import annotations

import importlib.util
import json
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("read_log", ROOT / "decomp" / "v4_recompute" / "read_log.py")
RL = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(RL)


def _lines(log):
    return [json.loads(x) for x in log.read_text(encoding="utf-8").splitlines()]


def test_the_line_is_on_disk_before_the_file_is_read(tmp_path, monkeypatch):
    log, f = tmp_path / "log.jsonl", tmp_path / "table.json"
    f.write_text("{}", encoding="utf-8")
    seen = []
    real = pathlib.Path.read_bytes

    def spy(self):
        seen.append(([x["단계"] for x in _lines(log)], self.name))
        return real(self)

    monkeypatch.setattr(pathlib.Path, "read_bytes", spy)
    RL.ReadLog(log, "test").before_read([f])
    assert seen == [(["before read"], "table.json")]
    assert [x["단계"] for x in _lines(log)] == ["before read", "after read"]


def test_a_log_that_cannot_be_written_stops_the_run(tmp_path):
    (tmp_path / "file").write_text("x", encoding="utf-8")
    with pytest.raises(RL.ReadLogError):
        RL.ReadLog(tmp_path / "file" / "log.jsonl", "test")


def test_a_log_lost_during_the_run_refuses_the_read(tmp_path, monkeypatch):
    log, f = tmp_path / "log.jsonl", tmp_path / "table.json"
    f.write_text("{}", encoding="utf-8")
    reads = []
    real = pathlib.Path.read_bytes
    monkeypatch.setattr(pathlib.Path, "read_bytes", lambda self: (reads.append(self.name), real(self))[1])
    r = RL.ReadLog(log, "test")
    r.path = tmp_path / "gone" / "log.jsonl"
    with pytest.raises(RL.ReadLogError):
        r.before_read([f])
    assert reads == []


def test_failure_control_a_read_without_the_log_leaves_no_line(tmp_path):
    log, f = tmp_path / "log.jsonl", tmp_path / "table.json"
    f.write_text("{}", encoding="utf-8")
    RL.ReadLog(log, "test")
    f.read_bytes()
    assert _lines(log) == []


def test_the_second_line_holds_the_sha256(tmp_path):
    import hashlib
    log, f = tmp_path / "log.jsonl", tmp_path / "table.json"
    f.write_bytes(b"abc")
    RL.ReadLog(log, "test").before_read([f, f])
    lines = _lines(log)
    assert len(lines) == 2 and lines[0]["바이트"] == 3 and lines[1]["sha256"] == hashlib.sha256(b"abc").hexdigest()


@pytest.mark.parametrize("script,args", [("run_subset.py", ["--pairs", "A~B", "--out", "OUT"]), ("compare.py", ["--selftest"])])
def test_the_scripts_stop_before_reading_when_the_log_cannot_be_written(tmp_path, script, args):
    (tmp_path / "file").write_text("x", encoding="utf-8")
    args = [str(tmp_path / "out") if a == "OUT" else a for a in args]
    r = subprocess.run([sys.executable, "-B", str(ROOT / "decomp" / "v4_recompute" / script), *args, "--read-log", str(tmp_path / "file" / "log.jsonl")],
                       capture_output=True, text=True, timeout=60)
    assert r.returncode != 0 and "cannot be written" in r.stderr
    assert not (tmp_path / "out" / "units").exists()
