#!/usr/bin/env python3

import os
import shlex
import subprocess
from pathlib import Path


CONFIG = {
    "launcher": "mpirun",
    "launcher_arguments": "-np 16",
    "executable": "../build_iwrap/reflectomITER",

    "input_template": "../inputs/input_IMAS.txt",
    "generated_input": "../inputs/input_IMAS_gen.txt",

    "working_directory": ".",
    "log_file": "reflectomITER_iwrap_python.log",
    "dry_run": False,

    "overrides": {
        "init.output_dir": "/home/pag/run_iwrap_test",
        "init.max_step": "3000",
        "init.final_time": "15e-9",
        "init.cfl": "0.45",
        "init.number_of_outputs": "5",
        "source.1.frequency": "3.7e9",
        "source.1.amplitude": "1e-3",
    },
}


STATE = {
    "step": 0,
    "status": "created",
    "last_command": "",
    "last_return_code": None,
}


def replace_or_append_parameter(lines: list[str], key: str, value: str) -> bool:
    for i, line in enumerate(lines):
        stripped = line.strip()

        if not stripped or stripped.startswith("#") or "=" not in line:
            continue

        left, _right = line.split("=", 1)

        if left.strip() == key:
            newline = "\n" if line.endswith("\n") else ""
            lines[i] = f"{key} = {value}{newline}"
            return True

    return False


def generate_input_file() -> str:
    template_path = Path(CONFIG["input_template"])
    generated_path = Path(CONFIG["generated_input"])

    if not template_path.exists():
        raise FileNotFoundError(f"Input template not found: {template_path}")

    generated_path.parent.mkdir(parents=True, exist_ok=True)

    lines = template_path.read_text(encoding="utf-8").splitlines(keepends=True)

    replaced = []
    appended = []

    for key, value in CONFIG["overrides"].items():
        value = str(value)

        if replace_or_append_parameter(lines, key, value):
            replaced.append((key, value))
        else:
            appended.append((key, value))

    if appended:
        if lines and not lines[-1].endswith("\n"):
            lines[-1] += "\n"

        lines.append("\n# ------------------------------------------------------------\n")
        lines.append("# iWrap-generated parameters not found in template\n")
        lines.append("# ------------------------------------------------------------\n")

        for key, value in appended:
            lines.append(f"{key} = {value}\n")

    generated_path.write_text("".join(lines), encoding="utf-8")

    report_path = generated_path.with_suffix(generated_path.suffix + ".report")

    with report_path.open("w", encoding="utf-8") as report:
        report.write("reflectomITER iWrap input-generation report\n")
        report.write(f"template: {template_path}\n")
        report.write(f"generated: {generated_path}\n\n")

        report.write("Replaced parameters:\n")
        for key, value in replaced:
            report.write(f"  {key} = {value}\n")

        report.write("\nAppended parameters:\n")
        for key, value in appended:
            report.write(f"  {key} = {value}\n")

    return str(generated_path)


def build_command() -> list[str]:
    input_file = generate_input_file()

    command = []

    if CONFIG["launcher"]:
        command.append(CONFIG["launcher"])

    if CONFIG["launcher_arguments"]:
        command.extend(shlex.split(CONFIG["launcher_arguments"]))

    command.extend([
        CONFIG["executable"],
        input_file,
    ])

    return command


def init_code():
    STATE["step"] = 0
    STATE["status"] = "initialized"
    STATE["last_command"] = ""
    STATE["last_return_code"] = None
    return 0, "reflectomITER Python actor initialized"


def clean_up():
    STATE["status"] = "finalized"
    return 0, "reflectomITER Python actor finalized"


def code_step():
    workdir = Path(CONFIG["working_directory"])
    log_file = workdir / CONFIG["log_file"]

    command = build_command()
    printable_command = " ".join(shlex.quote(x) for x in command)

    STATE["last_command"] = printable_command

    with log_file.open("a", encoding="utf-8") as log:
        log.write("\n[reflectomITER Python actor]\n")
        log.write(f"cwd     = {workdir}\n")
        log.write(f"command = {printable_command}\n")
        log.write(f"dry_run = {CONFIG['dry_run']}\n")

    if CONFIG["dry_run"]:
        STATE["step"] += 1
        STATE["status"] = "completed"
        STATE["last_return_code"] = 0
        return 0, f"DRY RUN: cd {workdir} && {printable_command}"

    try:
        with log_file.open("a", encoding="utf-8") as log:
            result = subprocess.run(
                command,
                cwd=str(workdir),
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
                text=True,
            )

        STATE["last_return_code"] = result.returncode

        if result.returncode == 0:
            STATE["step"] += 1
            STATE["status"] = "completed"
            return 0, "reflectomITER completed successfully"

        STATE["status"] = "failed"
        return result.returncode, f"reflectomITER failed with return code {result.returncode}"

    except FileNotFoundError as exc:
        STATE["status"] = "failed"
        STATE["last_return_code"] = 127
        return 127, f"launch failed: {exc}"


def get_code_state():
    state = (
        f"status={STATE['status']};"
        f"step={STATE['step']};"
        f"last_return_code={STATE['last_return_code']};"
        f"last_command={STATE['last_command']}"
    )
    return state, 0, "state returned"


def restore_code_state(state: str):
    STATE["status"] = state
    return 0, "state restored"


def get_timestamp_cpp():
    return float(STATE["step"]), 0, "timestamp returned"


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