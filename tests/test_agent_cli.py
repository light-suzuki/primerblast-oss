import argparse
import io
import json
import subprocess
import sys
import pytest

from primerblast_oss import agent_cli, workflows
from primerblast_oss.webapp import server


def test_schema_is_discoverable_without_gui():
    process = subprocess.run([sys.executable, "-m", "primerblast_oss", "agent", "schema"],
                             capture_output=True, encoding="utf-8")
    assert process.returncode == 0
    schema = json.loads(process.stdout)
    assert schema["starts_gui"] is False
    assert {operation["name"] for operation in schema["operations"]} == set(server.HANDLERS)


def test_agent_and_gui_share_execution_and_keep_stdout_json(monkeypatch, capsys):
    calls = []
    def handler(params):
        calls.append(params)
        print("engine diagnostic")
        return {"mode": "check", "search_status": "incomplete", "results": []}
    monkeypatch.setitem(workflows.HANDLERS, "check", handler)
    monkeypatch.setattr(sys, "stdin", io.StringIO(json.dumps({"operation": "check", "params": {"db": ["synthetic"]}, "request_id": "r1"})))
    assert agent_cli.command(argparse.Namespace(action="run", input="-", allow_db_write=False)) == 0
    output = capsys.readouterr()
    result = json.loads(output.out)
    assert "engine diagnostic" in output.err
    assert result["request_id"] == "r1"
    assert result["result"]["search_status"] == "incomplete"
    jobs = server.JobManager()
    monkeypatch.setattr(server.threading.Thread, "start", lambda self: self.run())
    job = jobs.get(jobs.submit("check", {"db": ["synthetic"]}))
    assert job["status"] == "done" and job["result"] == result["result"]
    assert len(calls) == 2


def test_database_write_is_blocked_before_engine_execution(monkeypatch, capsys):
    monkeypatch.setitem(workflows.HANDLERS, "makedb", lambda _: (_ for _ in ()).throw(AssertionError("must not execute")))
    monkeypatch.setattr(sys, "stdin", io.StringIO('{"operation":"makedb","params":{}}'))
    assert agent_cli.command(argparse.Namespace(action="run", input="-", allow_db_write=False)) == 2
    assert "allow-db-write" in json.loads(capsys.readouterr().out)["error"]["message"]


def test_bad_request_returns_json_and_nonzero_exit():
    process = subprocess.run([sys.executable, "-m", "primerblast_oss", "agent", "run"],
                             input='{"operation":"unknown","params":{}}', capture_output=True, encoding="utf-8")
    assert process.returncode == 2 and json.loads(process.stdout)["ok"] is False


def test_application_and_agent_do_not_import_gui_server():
    result = subprocess.run([sys.executable, "-c", "import primerblast_oss.agent_cli, sys; assert 'primerblast_oss.webapp.server' not in sys.modules"],
                            capture_output=True, encoding="utf-8")
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("result", [float("nan"), {"not_json": {1, 2}}])
def test_unencodable_engine_result_returns_json_error(monkeypatch, capsys, result):
    monkeypatch.setitem(workflows.HANDLERS, "check", lambda _: result)
    monkeypatch.setattr(sys, "stdin", io.StringIO('\ufeff{"operation":"check","params":{},"request_id":"r1"}'))
    assert agent_cli.command(argparse.Namespace(action="run", input="-", allow_db_write=False)) == 1
    output = json.loads(capsys.readouterr().out)
    assert output["error"]["code"] == "execution_failed" and output["request_id"] == "r1"


@pytest.mark.parametrize("number", ["NaN", "Infinity", "-Infinity"])
def test_nonfinite_input_is_rejected_before_execution(monkeypatch, capsys, number):
    monkeypatch.setitem(workflows.HANDLERS, "check", lambda _: pytest.fail("must not execute"))
    monkeypatch.setattr(sys, "stdin", io.StringIO('{"operation":"check","params":{"value":' + number + '}}'))
    assert agent_cli.command(argparse.Namespace(action="run", input="-", allow_db_write=False)) == 2
    assert json.loads(capsys.readouterr().out)["error"]["code"] == "invalid_request"
