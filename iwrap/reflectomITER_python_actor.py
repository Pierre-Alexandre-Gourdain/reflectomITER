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

def _env_optional(name: str) -> str | None:
    value = os.environ.get(name)
    if value is None:
        return None
    value = value.strip()
    return value if value else None
    
def _replace_or_append_parameter(lines: list[str], key: str, value: str) -> bool:
    """
    Replace a ParmParse-style assignment:

        key = old_value

    while preserving comments and unrelated lines.

    Returns True if replacement happened, False otherwise.
    """
    prefix = f"{key}"

    for i, line in enumerate(lines):
        stripped = line.strip()

        if not stripped or stripped.startswith("#"):
            continue

        if "=" not in line:
            continue

        left, _right = line.split("=", 1)

        if left.strip() == prefix:
            newline = "\n" if line.endswith("\n") else ""
            lines[i] = f"{key} = {value}{newline}"
            return True

    return False


def generate_input_file() -> str:
    template_path = Path(
        os.environ.get(
            "REFLECTOMITER_INPUT_TEMPLATE",
            "templates/input.reflectomITER.template",
        )
    )

    output_input_path = Path(
        os.environ.get(
            "REFLECTOMITER_GENERATED_INPUT",
            "generated/inputs.iwrap_reflectomITER",
        )
    )

    if not template_path.exists():
        raise FileNotFoundError(f"Input template not found: {template_path}")

    output_input_path.parent.mkdir(parents=True, exist_ok=True)

    lines = template_path.read_text(encoding="utf-8").splitlines(keepends=True)

    source_index = os.environ.get("REFLECTOMITER_SOURCE_INDEX", "1").strip()

    mapping = {
        "REFLECTOMITER_OUTPUT_DIR": "init.output_dir",
        "REFLECTOMITER_MAX_STEP": "init.max_step",
        "REFLECTOMITER_FINAL_TIME": "init.final_time",
        "REFLECTOMITER_CFL": "init.cfl",
        "REFLECTOMITER_NUMBER_OF_OUTPUTS": "init.number_of_outputs",
        "REFLECTOMITER_SOURCE_FREQUENCY": f"source.{source_index}.frequency",
        "REFLECTOMITER_SOURCE_AMPLITUDE": f"source.{source_index}.amplitude",
    }

    replaced: list[tuple[str, str]] = []
    appended: list[tuple[str, str]] = []

    for env_name, input_key in mapping.items():
        value = _env_optional(env_name)
        if value is None:
            continue

        did_replace = _replace_or_append_parameter(lines, input_key, value)

        if did_replace:
            replaced.append((input_key, value))
        else:
            appended.append((input_key, value))

    if appended:
        if lines and not lines[-1].endswith("\n"):
            lines[-1] += "\n"

        lines.append("\n# ------------------------------------------------------------\n")
        lines.append("# iWrap-generated parameters not found in template\n")
        lines.append("# ------------------------------------------------------------\n")

        for key, value in appended:
            lines.append(f"{key} = {value}\n")

    output_input_path.write_text("".join(lines), encoding="utf-8")

    # Write a small sidecar report for debugging.
    report_path = output_input_path.with_suffix(output_input_path.suffix + ".report")
    with report_path.open("w", encoding="utf-8") as report:
        report.write("reflectomITER iWrap input-generation report\n")
        report.write(f"template: {template_path}\n")
        report.write(f"generated: {output_input_path}\n\n")

        report.write("Replaced parameters:\n")
        for key, value in replaced:
            report.write(f"  {key} = {value}\n")

        report.write("\nAppended parameters:\n")
        for key, value in appended:
            report.write(f"  {key} = {value}\n")

    return str(output_input_path)   

   
def _command() -> list[str]:
    executable = os.environ.get("REFLECTOMITER_EXECUTABLE", "../build/reflectomITER")

    if _env_bool("REFLECTOMITER_GENERATE_INPUT"):
        input_file = generate_input_file()
    else:
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