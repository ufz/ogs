import sys
from pathlib import Path

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: compare_binaries.py <cmake list>")
        sys.exit(1)

    cmake_binaries = sys.argv[1].split(";")

    script_dir = Path(__file__).resolve().parent
    parent_dir = (script_dir / "../../Applications/Python").resolve()
    sys.path.insert(0, str(parent_dir))

    from ogs._internal.binaries_list import binaries_list

    cmake_set = set(cmake_binaries)
    python_set = set(binaries_list)
    missing_in_python = cmake_set - python_set

    if missing_in_python:
        print("Missing in Pythons binaries_list:", sorted(missing_in_python))
        sys.exit(1)

    sys.exit(0)
