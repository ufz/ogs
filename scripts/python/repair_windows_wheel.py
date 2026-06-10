import argparse
import os
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--wheel-dir", required=True)
    parser.add_argument("wheel")
    return parser.parse_args()


def main():
    args = parse_args()
    wheel = Path(args.wheel).resolve()
    wheel_dir = Path(args.wheel_dir).resolve()

    with tempfile.TemporaryDirectory() as tmp:
        extract_dir = Path(tmp) / "wheel"
        with zipfile.ZipFile(wheel) as whl:
            whl.extractall(extract_dir)

        add_paths = [
            str(path)
            for path in (extract_dir / "ogs", extract_dir / "bin")
            if path.is_dir()
        ]

        command = [
            sys.executable,
            "-m",
            "delvewheel",
            "repair",
            "--add-path",
            os.pathsep.join(add_paths),
            "-w",
            str(wheel_dir),
            "-v",
            str(wheel),
        ]
        subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
