from __future__ import annotations

import sys
import os

sys.path.append(os.path.join("Applications", "Python"))
from ogs._internal.provide_ogs_cli_tools_via_wheel import pyproject_get_scripts

__all__ = ["dynamic_metadata"]


def __dir__() -> list[str]:
    return __all__


def dynamic_metadata(
    fields: frozenset[str],
    settings: dict[str, object] | None = None,
) -> dict[str, str | dict[str, str | None]]:
    if fields != {"scripts"}:
        msg = "Only the 'scripts' field is supported"
        raise ValueError(msg)

    if settings:
        msg = "No inline configuration is supported"
        raise ValueError(msg)

    # Structure for entry-points
    # eps = {
    #    "console_scripts": {"ogs": "ogs._internal.provide_ogs_cli_tools_via_wheel:ogs"}
    # }
    # scripts = {"ogs": "ogs._internal.provide_ogs_cli_tools_via_wheel:ogs"}

    return {"scripts": pyproject_get_scripts()}
