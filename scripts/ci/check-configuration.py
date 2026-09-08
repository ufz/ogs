#!/usr/bin/env python3
"""Check the repository's literal CI and CMake preset references."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

LOCAL_INCLUDE = re.compile(r"^\s*-\s*local:\s*[\"']?([^\"'#\s]+)")
CI_PRESET = re.compile(r"^\s+(CMAKE_PRESET|CTEST_PRESET):\s*[\"']?([^\"'#\s]+)")

# Keep exceptional build/configure aliases explicit and reviewable. Dynamic
# references and empty defaults are intentionally outside this check.
DOCUMENTED_BUILD_ALIASES: dict[str, str] = {}


def ci_files(root: Path) -> list[Path]:
    return [
        root / ".gitlab-ci.yml",
        *sorted((root / "scripts/ci").rglob("*.yml")),
    ]


def check_local_includes(root: Path) -> list[str]:
    errors = []
    for filename in ci_files(root):
        for line_number, line in enumerate(filename.read_text().splitlines(), 1):
            match = LOCAL_INCLUDE.match(line)
            if not match:
                continue
            include = match.group(1).lstrip("/")
            if any(character in include for character in "*?["):
                matches = list(root.glob(include))
                if not matches:
                    errors.append(
                        f"{filename}:{line_number}: include glob matches no files: {include}"
                    )
            elif not (root / include).is_file():
                errors.append(
                    f"{filename}:{line_number}: missing local include: {include}"
                )
    return errors


def check_presets(root: Path) -> list[str]:
    filename = root / "CMakePresets.json"
    data = json.loads(filename.read_text())
    errors = []
    configure = {preset["name"] for preset in data.get("configurePresets", [])}
    builds = {preset["name"]: preset for preset in data.get("buildPresets", [])}
    tests = {preset["name"]: preset for preset in data.get("testPresets", [])}

    for preset in builds.values():
        reference = preset.get("configurePreset")
        if reference and reference not in configure:
            errors.append(
                f"{filename}: build preset {preset['name']} references missing "
                f"configure preset {reference}"
            )
        if (
            preset["name"] in configure
            and reference != preset["name"]
            and DOCUMENTED_BUILD_ALIASES.get(preset["name"]) != reference
        ):
            errors.append(
                f"{filename}: build preset {preset['name']} maps to {reference} "
                "instead of its same-named configure preset"
            )

    for preset in tests.values():
        reference = preset.get("configurePreset")
        if reference and reference not in configure:
            errors.append(
                f"{filename}: test preset {preset['name']} references missing "
                f"configure preset {reference}"
            )
    for filename in ci_files(root):
        for line_number, line in enumerate(filename.read_text().splitlines(), 1):
            match = CI_PRESET.match(line)
            if not match or match.group(2).startswith("$"):
                continue
            preset_type, name = match.groups()
            valid_names = configure if preset_type == "CMAKE_PRESET" else tests
            if name and name not in valid_names:
                errors.append(
                    f"{filename}:{line_number}: {preset_type} references missing "
                    f"preset {name}"
                )
    return errors


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    errors = [*check_local_includes(root), *check_presets(root)]
    if errors:
        print("CI/CMake configuration check failed:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("CI/CMake configuration references are consistent.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
