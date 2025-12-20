# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import os
import subprocess
import sys
from pathlib import Path

root_dir = Path(__file__).resolve().parent.parent.parent
os.chdir(root_dir / "web" / "content")
notebooks = Path().glob("**/*.ipynb")
exit_code = 0

for notebook in notebooks:
    nb = Path(notebook)
    if nb.parent.stem != nb.stem:
        print(
            f"Error for {notebook}: Notebook must have the same name as its parent folder!"
        )
        exit_code = 1
        continue
    template = root_dir / "Tests/Data/Notebooks/nbconvert_templates/collapsed.md.j2"
    subprocess.run(
        [
            "jupyter",
            "nbconvert",
            "--to",
            "markdown",
            f"--template-file={template}",
            "--output=index",
            notebook,
        ],
        check=True,
    )

sys.exit(exit_code)
