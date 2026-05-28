#!/usr/bin/env python3

import os
import shlex
import subprocess
from pathlib import Path

_state = {
    "step": 0,
    "status": "created",
    "last_command": "",
    "last_return_code": None,
}


def _env_bool(name: str) -> bool:
    return os.environ.get(name, "").strip().lower() in {"1", "true", "yes", "on"}


def _command() -> list[str]:
    executable = os.environ.get("REFLECTOMITER_EXECUTABLE", "../build/reflectomITER")
    input_file = os.environ.get("REFLECTOMITER_INPUT", "inputs.iwrap_reflectomITER")
    launcher = os.environ.get("REFLECTOMITER_LAUNCHER", "")
    launcher_args = os.environ.get("REFLECTOMITER_LAUNCHER_ARGS", "")

    cmd: list[str] = []

    if launcher:
        cmd.append(launcher)

    if launcher_args:
        cmd.extend(shlex.split(launcher_args))

    cmd.extend([executable, input_file])
    return cmd


def init_code():
    _state["step"] = 0
    _state["status"] = "initialized"
    _state["last_command"] = ""
    _state["last_return_code"] = None
    return 0, "reflectomITER Python actor initialized"


def clean_up():
    _state["status"] = "finalized"
    return 0, "reflectomITER Python actor finalized"


def code_step():
    workdir = Path(os.environ.get("REFLECTOMITER_WORKDIR", "."))
    log_file = workdir / os.environ.get("REFLECTOMITER_IWRAP_LOG", "reflectomITER_iwrap_python.log")
    dry_run = _env_bool("REFLECTOMITER_DRY_RUN")

    cmd = _command()
    printable = " ".join(shlex.quote(x) for x in cmd)
    _state["last_command"] = printable

    with log_file.open("a", encoding="utf-8") as log:
        log.write("\n[reflectomITER Python actor]\n")
        log.write(f"cwd     = {workdir}\n")
        log.write(f"command = {printable}\n")
        log.write(f"dry_run = {dry_run}\n")

    if dry_run:
        _state["step"] += 1
        _state["status"] = "completed"
        _state["last_return_code"] = 0
        return 0, f"DRY RUN: cd {workdir} && {printable}"

    try:
        with log_file.open("a", encoding="utf-8") as log:
            result = subprocess.run(
                cmd,
                cwd=str(workdir),
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
                text=True,
            )

        _state["last_return_code"] = result.returncode

        if result.returncode == 0:
            _state["step"] += 1
            _state["status"] = "completed"
            return 0, "reflectomITER completed successfully"

        _state["status"] = "failed"
        return result.returncode, f"reflectomITER failed with return code {result.returncode}"

    except FileNotFoundError as exc:
        _state["status"] = "failed"
        _state["last_return_code"] = 127
        return 127, f"launch failed: {exc}"


def get_code_state():
    state = (
        f"status={_state['status']};"
        f"step={_state['step']};"
        f"last_return_code={_state['last_return_code']};"
        f"last_command={_state['last_command']}"
    )
    return state, 0, "state returned"


def restore_code_state(state: str):
    _state["status"] = state
    return 0, "state restored"


def get_timestamp_cpp():
    return float(_state["step"]), 0, "timestamp returned"


def main() -> int:
    status, message = init_code()
    print(f"init: {status} | {message}")

    status, message = code_step()
    print(f"step: {status} | {message}")

    state, status, message = get_code_state()
    print(f"state: {state}")

    status, message = clean_up()
    print(f"cleanup: {status} | {message}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())