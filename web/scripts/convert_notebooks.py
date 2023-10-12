import os
import subprocess
import sys
from pathlib import Path

content_dir = Path(__file__).resolve().parent.parent / "content"
os.chdir(content_dir)
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
    template = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "../../Tests/Data/Notebooks/nbconvert_templates/collapsed.md.j2",
    )
    subprocess.run(
        [
            "jupyter",
            "nbconvert",
            "--to",
            "markdown",
            f"--template-file={template}",
            "--output=index",
            notebook,
        ]
    )

sys.exit(exit_code)
