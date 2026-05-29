#!/usr/bin/env python3
"""
Convert Doxygen-style MathJax documentation to GitHub Markdown math.

Conversions:

    \\f$ ... \\f$           -> $`...`$
    \\f[ ... \\f]           -> ```math ... ```
    \\f{align}{ ... \\f}    -> ```math \\begin{aligned} ... \\end{aligned} ```
    \\f{equation}{ ... \\f} -> ```math ... ```

Edit INPUT_FILE and OUTPUT_FILE below.
"""

from pathlib import Path
import re


INPUT_FILE = Path("../src/main_page.md")
OUTPUT_FILE = Path("../README.md")


def convert_inline_math(text: str) -> str:
    """
    Convert Doxygen inline math:

        \\f$ E = mc^2 \\f$

    to GitHub's robust inline math syntax:

        $`E = mc^2`$
    """

    pattern = re.compile(
        r"\\f\$(.*?)\\f\$",
        flags=re.DOTALL,
    )

    def repl(match: re.Match) -> str:
        content = match.group(1).strip()
        return f"$`{content}`$"

    return pattern.sub(repl, text)


def convert_display_bracket_math(text: str) -> str:
    """
    Convert Doxygen display math:

        \\f[
        E = mc^2
        \\f]

    to GitHub fenced math:

        ```math
        E = mc^2
        ```
    """

    pattern = re.compile(
        r"\\f\[(.*?)\\f\]",
        flags=re.DOTALL,
    )

    def repl(match: re.Match) -> str:
        content = match.group(1).strip()
        return f"\n```math\n{content}\n```\n"

    return pattern.sub(repl, text)


def convert_doxygen_environment_math(text: str) -> str:
    """
    Convert Doxygen named math environments.

    Examples:

        \\f{align}{
        a &= b \\\\
        c &= d
        \\f}

    becomes:

        ```math
        \\begin{aligned}
        a &= b \\\\
        c &= d
        \\end{aligned}
        ```
    """

    pattern = re.compile(
        r"\\f\{([A-Za-z*]+)\}\{(.*?)\\f\}",
        flags=re.DOTALL,
    )

    def repl(match: re.Match) -> str:
        env = match.group(1).strip()
        content = match.group(2).strip()

        if env in {"align", "align*"}:
            body = (
                "\\begin{aligned}\n"
                f"{content}\n"
                "\\end{aligned}"
            )
        elif env in {"equation", "equation*"}:
            body = content
        else:
            body = (
                f"\\begin{{{env}}}\n"
                f"{content}\n"
                f"\\end{{{env}}}"
            )

        return f"\n```math\n{body}\n```\n"

    return pattern.sub(repl, text)


def normalize_blank_lines(text: str) -> str:
    """
    Avoid excessive blank lines introduced by block replacements.
    """

    return re.sub(r"\n{4,}", "\n\n\n", text)


def convert_doxygen_to_github(text: str) -> str:
    """
    Apply all conversions.

    Order matters:
    1. Named environments
    2. Display bracket math
    3. Inline math
    """

    text = convert_doxygen_environment_math(text)
    text = convert_display_bracket_math(text)
    text = convert_inline_math(text)
    text = normalize_blank_lines(text)
    return text


def main() -> None:
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

    source = INPUT_FILE.read_text(encoding="utf-8")
    converted = convert_doxygen_to_github(source)

    OUTPUT_FILE.write_text(converted, encoding="utf-8")

    print("Converted:")
    print(f"  input : {INPUT_FILE}")
    print(f"  output: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()