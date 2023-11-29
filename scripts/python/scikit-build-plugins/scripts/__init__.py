from __future__ import annotations

import sys
from pathlib import Path

sys.path.append(str(Path("Applications").joinpath("Python").absolute()))
# ruff: noqa: E402
from ogs._internal.provide_ogs_cli_tools_via_wheel import (
    pyproject_get_binaries,
)

__all__ = ["dynamic_metadata"]


def __dir__() -> list[str]:
    return __all__


def dynamic_metadata(
    field: str,
    settings: dict[str, object] | None = None,
) -> str:
    if field != "scripts":
        msg = "Only the 'scripts' field is supported"
        raise ValueError(msg)

    if settings:
        msg = "No inline configuration is supported"
        raise ValueError(msg)

    return pyproject_get_binaries()
