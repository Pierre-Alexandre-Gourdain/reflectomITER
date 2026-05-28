def make_source(index: int, base_source: dict[str, str], updates: dict[str, str]) -> dict[str, str]:
    source = dict(base_source)
    source.update(updates)
    source["index"] = str(index)
    return source


def source_to_overrides(source: dict[str, str]) -> dict[str, str]:
    index = source["index"]
    prefix = f"source.{index}"

    overrides = {}

    for key, value in source.items():
        if key == "index":
            continue

        overrides[f"{prefix}.{key}"] = str(value)

    return overrides


def sources_to_overrides(sources: list[dict[str, str]]) -> dict[str, str]:
    overrides = {
        "init.nsources": str(len(sources)),
    }

    for source in sources:
        overrides.update(source_to_overrides(source))

    return overrides