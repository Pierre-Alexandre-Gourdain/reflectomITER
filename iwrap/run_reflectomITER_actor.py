#!/usr/bin/env python3

from pathlib import Path

from reflectomiter_actor.actor import ReflectomITERActor
from reflectomiter_actor.config import ActorConfig


def main() -> int:
    config = ActorConfig(
        launcher="mpirun",
        launcher_arguments="-np 1",
        executable="../build/reflectomITER",

        input_template=Path("../inputs/input_IMAS.txt"),
        generated_input=Path("../inputs/input_IMAS_gen.txt"),

        working_directory=Path("."),
        log_file=Path("reflectomITER_iwrap_python.log"),
        dry_run=False,

        overrides={
            "init.output_dir": "/home/pag/run_iwrap_test",
            "init.max_step": "100",
            "init.final_time": "15e-9",
            "init.cfl": "0.45",
            "init.number_of_outputs": "5",
            "source.1.frequency": "3.7e9",
            "source.1.amplitude": "1e-3",
        },
    )

    actor = ReflectomITERActor(config)

    status, message = actor.init_code()
    print(f"init: {status} | {message}")

    status, message = actor.code_step()
    print(f"step: {status} | {message}")

    state, status_state, message_state = actor.get_code_state()
    print(f"state: {state}")

    status_cleanup, message_cleanup = actor.clean_up()
    print(f"cleanup: {status_cleanup} | {message_cleanup}")

    return status


if __name__ == "__main__":
    raise SystemExit(main())