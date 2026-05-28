from .config import ActorConfig
from .input_generation import generate_input_file
from .launcher import build_command, run_command


class ReflectomITERActor:
    def __init__(self, config: ActorConfig):
        self.config = config
        self.step = 0
        self.status = "created"
        self.last_command = ""
        self.last_return_code = None

    def init_code(self) -> tuple[int, str]:
        self.step = 0
        self.status = "initialized"
        self.last_command = ""
        self.last_return_code = None
        return 0, "reflectomITER Python actor initialized"

    def clean_up(self) -> tuple[int, str]:
        self.status = "finalized"
        return 0, "reflectomITER Python actor finalized"

    def code_step(self) -> tuple[int, str]:
        input_file = generate_input_file(
            template_path=self.config.input_template,
            generated_path=self.config.generated_input,
            overrides=self.config.overrides,
        )

        command = build_command(self.config, input_file)
        self.last_command = " ".join(command)

        status_code, status_message = run_command(self.config, command)
        self.last_return_code = status_code

        if status_code == 0:
            self.step += 1
            self.status = "completed"
        else:
            self.status = "failed"

        return status_code, status_message

    def get_code_state(self) -> tuple[str, int, str]:
        state = (
            f"status={self.status};"
            f"step={self.step};"
            f"last_return_code={self.last_return_code};"
            f"last_command={self.last_command}"
        )

        return state, 0, "state returned"

    def restore_code_state(self, state: str) -> tuple[int, str]:
        self.status = state
        return 0, "state restored"

    def get_timestamp_cpp(self) -> tuple[float, int, str]:
        return float(self.step), 0, "timestamp returned"