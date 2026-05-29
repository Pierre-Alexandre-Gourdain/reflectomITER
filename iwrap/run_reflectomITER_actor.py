#!/usr/bin/env python3

from pathlib import Path

from reflectomITER_actor.actor import ReflectomITERActor
from reflectomITER_actor.config import ActorConfig
from reflectomITER_actor.sources import make_source, sources_to_overrides


BASE_RUN_OVERRIDES = {
    "init.max_step": "-1",
    "init.final_time": "50e-9",
    "init.cfl": "0.45",
    "init.number_of_outputs": "50",
}


BASE_SOURCE = {
    "enabled": "1",
    "type": "J",
    "position": "8.5 0.0 0.5",
    "direction": "0.0 1.0 0.0",
    "shape": "0.05 1.0 1.0",
    "amplitude": "1e-3",
    "frequency": "37e9",
    "phase": "0.0",
    "width": ".05",
    "profile": "gaussian",
    "supergaussian_order": "4.0",
    "t_on": "1e-9",
    "t_off": "11e-9",
    "t_rise": ".3e-9",
}


SOURCES = [
    make_source(
        index=1,
        base_source=BASE_SOURCE,
        updates={
            "frequency": "36.5e9",
        },
    ),
    make_source(
        index=2,
        base_source=BASE_SOURCE,
        updates={
            "frequency": "37e9",
        },
    ),
    make_source(
        index=3,
        base_source=BASE_SOURCE,
        updates={
            "frequency": "38.5e9",
        },
    ),
    make_source(
        index=4,
        base_source=BASE_SOURCE,
        updates={
            "frequency": "40e9",
        },
    ),
    make_source(
        index=5,
        base_source=BASE_SOURCE,
        updates={
            "frequency": "41.5e9",
        },
    ),
]


def make_config() -> ActorConfig:
    overrides = {
        **BASE_RUN_OVERRIDES,
        **sources_to_overrides(SOURCES),
    }

    return ActorConfig(
        launcher="mpirun",
        launcher_arguments="-np 1",
        executable="../build/reflectomITER",

        input_template=Path("../inputs/input_IMAS.txt"),
        generated_input=Path("../inputs/input_IMAS_multisource_gen.txt"),

        working_directory=Path("."),
        log_file=Path("reflectomITER_iwrap_multisource.log"),
        dry_run=False,

        overrides=overrides,
    )


def main() -> int:
    actor = ReflectomITERActor(make_config())

    status, message = actor.init_code()
    print(f"init: {status} | {message}")

    status, message = actor.code_step()
    print(f"step: {status} | {message}")

    state, _, _ = actor.get_code_state()
    print(f"state: {state}")

    status, message = actor.clean_up()
    print(f"cleanup: {status} | {message}")

    return status


if __name__ == "__main__":
    raise SystemExit(main())