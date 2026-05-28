from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ActorConfig:
    launcher: str = "mpirun"
    launcher_arguments: str = "-np 1"
    executable: str = "../build/reflectomITER"

    input_template: Path = Path("../inputs/input_IMAS.txt")
    generated_input: Path = Path("../inputs/input_IMAS_gen.txt")

    working_directory: Path = Path(".")
    log_file: Path = Path("reflectomITER_iwrap_python.log")
    dry_run: bool = False

    overrides: dict[str, str] = field(default_factory=lambda: {
        "init.output_dir": "../run_iwrap_test",
        "init.max_step": "100",
        "init.final_time": "15e-9",
        "init.cfl": "0.45",
        "init.number_of_outputs": "5",
        "source.1.frequency": "3.7e9",
        "source.1.amplitude": "1e-3",
    })