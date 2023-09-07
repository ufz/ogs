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
    print(f"Converting {notebook} ...")
    subprocess.run(
        f"nb2hugo --site-dir .. --section {nb.parent.parent} {notebook}",
        shell=True,
        check=True,
    )

sys.exit(exit_code)
