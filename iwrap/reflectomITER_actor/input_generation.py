from pathlib import Path


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


def generate_input_file(
    template_path: Path,
    generated_path: Path,
    overrides: dict[str, str],
) -> Path:
    if not template_path.exists():
        raise FileNotFoundError(f"Input template not found: {template_path}")

    generated_path.parent.mkdir(parents=True, exist_ok=True)

    lines = template_path.read_text(encoding="utf-8").splitlines(keepends=True)

    replaced: list[tuple[str, str]] = []
    appended: list[tuple[str, str]] = []

    for key, value in overrides.items():
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

    write_report(
        report_path=generated_path.with_suffix(generated_path.suffix + ".report"),
        template_path=template_path,
        generated_path=generated_path,
        replaced=replaced,
        appended=appended,
    )

    return generated_path


def write_report(
    report_path: Path,
    template_path: Path,
    generated_path: Path,
    replaced: list[tuple[str, str]],
    appended: list[tuple[str, str]],
) -> None:
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