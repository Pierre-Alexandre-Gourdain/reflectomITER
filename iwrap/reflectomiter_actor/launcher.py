import shlex
import subprocess
from pathlib import Path

from .config import ActorConfig


def build_command(config: ActorConfig, input_file: Path) -> list[str]:
    command: list[str] = []

    if config.launcher:
        command.append(config.launcher)

    if config.launcher_arguments:
        command.extend(shlex.split(config.launcher_arguments))

    command.extend([
        config.executable,
        str(input_file),
    ])

    return command


def run_command(config: ActorConfig, command: list[str]) -> tuple[int, str]:
    workdir = config.working_directory
    log_path = workdir / config.log_file

    printable_command = " ".join(shlex.quote(x) for x in command)

    with log_path.open("a", encoding="utf-8") as log:
        log.write("\n[reflectomITER Python actor]\n")
        log.write(f"cwd     = {workdir}\n")
        log.write(f"command = {printable_command}\n")
        log.write(f"dry_run = {config.dry_run}\n")

    if config.dry_run:
        return 0, f"DRY RUN: cd {workdir} && {printable_command}"

    try:
        with log_path.open("a", encoding="utf-8") as log:
            result = subprocess.run(
                command,
                cwd=str(workdir),
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
                text=True,
            )

        if result.returncode == 0:
            return 0, "reflectomITER completed successfully"

        return result.returncode, f"reflectomITER failed with return code {result.returncode}"

    except FileNotFoundError as exc:
        return 127, f"launch failed: {exc}"