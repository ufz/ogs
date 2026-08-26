# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import argparse
import os
import shutil
import subprocess
import sys
from datetime import timedelta
from pathlib import Path
from timeit import default_timer as timer

import jupytext
import nbformat
import papermill
import toml
from nbclient.exceptions import DeadKernelError
from nbconvert import HTMLExporter
from nbconvert.preprocessors import CellExecutionError


def coverage_enabled():
    value = os.getenv("OGS_COVERAGE_PYTHON", "").lower()
    return value in {"1", "true", "yes", "on"}


def petsc_enabled():
    value = os.getenv("OGS_USE_PETSC", "").lower()
    return value in {"1", "true", "yes", "on"}


def get_website_output_path(notebook, exec_notebook_file):
    if "Tests/Data" not in str(exec_notebook_file):
        return None

    first_cell = notebook.cells[0]
    lines = first_cell.source.splitlines()
    toml_begin = lines.index("+++")
    toml_end = max(loc for loc, val in enumerate(lines) if val == "+++")
    toml_lines = lines[toml_begin + 1 : toml_end]
    parsed_frontmatter = toml.loads("\n".join(toml_lines))
    return (
        Path(build_dir)
        / Path("web/content/docs/benchmarks")
        / Path(parsed_frontmatter["web_subsection"])
        / exec_notebook_file.stem.lower()
    )


def save_to_website(exec_notebook_file):
    output_path_arg = ""
    notebook = nbformat.read(exec_notebook_file, as_version=4)
    output_path = get_website_output_path(notebook, exec_notebook_file)
    if output_path:
        output_path_arg = f"--output-dir={Path(output_path)}"

    template = Path(__file__).parent.resolve() / "nbconvert_templates/collapsed.md.j2"
    subprocess.run(
        [
            "jupyter",
            "nbconvert",
            "--to",
            "markdown",
            f"--template-file={template}",
            "--output=index",
            output_path_arg,
            exec_notebook_file,
        ],
        check=True,
    )

    if not output_path:
        return

    Path(output_path).mkdir(parents=True, exist_ok=True)

    to_copy_path = Path(output_path) / "."

    if is_jupytext:
        shutil.copy(exec_notebook_file, to_copy_path)
        print(f"Copying ${exec_notebook_file} to {to_copy_path}")

    for subfolder in ["figures", "images"]:
        figures_path = (Path(notebook_file_path).parent / subfolder).resolve()
        symlink_figures_path = to_copy_path / subfolder
        if figures_path.exists():
            print(
                f"{subfolder} folder detected, copying {figures_path} to "
                f"{symlink_figures_path}"
            )
            shutil.copytree(figures_path, symlink_figures_path, dirs_exist_ok=True)


def check_and_modify_frontmatter():
    success = True
    # Check frontmatter has its own cell
    first_cell = nb["cells"][0]
    if (
        first_cell.cell_type == "markdown"
        and first_cell.source.startswith("+++")
        and not first_cell.source.endswith("+++")
    ):
        print(
            f"Error: {notebook_filename} notebook metadata is not a separate cell (in markdown: separate by two newlines)!"
        )
        success = False

    # Modify metadata
    first_cell = nb["cells"][0]
    if first_cell.source.startswith("---"):
        print(
            f"Error: {notebook_filename} frontmatter is not in TOML format! Use +++ delimitiers!"
        )
        success = False
    # Add metadata used by the Hugo layout to render the notebook banner.
    repo = "https://gitlab.opengeosys.org/ogs/ogs"
    branch = "master"
    binder_tag = "6.5.8-0.8.2"
    petsc_prefix = ""
    if petsc_enabled():
        petsc_prefix = "petsc-"
    if "CI_MERGE_REQUEST_SOURCE_PROJECT_URL" in os.environ:
        repo = os.environ["CI_MERGE_REQUEST_SOURCE_PROJECT_URL"]
        branch = os.environ["CI_MERGE_REQUEST_SOURCE_BRANCH_NAME"]
    binder_link = f"https://binder.opengeosys.org/v2/gh/bilke/binder-ogs-requirements/{petsc_prefix}{binder_tag}?urlpath=git-pull%3Frepo={repo}%26urlpath=lab/tree/ogs/{notebook_file_path_relative}%26branch={branch}%26depth=1"
    metadata = (
        "notebook = true\n"
        f'notebook_source_url = "{repo}/-/blob/{branch}/{notebook_file_path_relative}"\n'
        f'notebook_source_name = "{notebook_filename}"\n'
        f'notebook_binder_url = "{binder_link}"\n'
    )
    if is_jupytext:
        metadata += f'notebook_download_file = "{Path(convert_notebook_file).name}"\n'
    first_cell.source = first_cell.source.replace("+++\n", f"+++\n{metadata}", 1)

    return success


# Script arguments
parser = argparse.ArgumentParser(description="Jupyter notebook testrunner.")
parser.add_argument("notebooks", metavar="N", nargs="+", help="Notebooks to test.")
parser.add_argument("--out", default="./", help="Output directory.")
parser.add_argument(
    "--hugo", action="store_true", help="Convert successful notebooks to web site."
)
parser.add_argument("--hugo-out", default="web", help="Hugo output directory.")
args = parser.parse_args()

# Path setup
testrunner_script_path = Path(__file__).parent.resolve()
ogs_source_path = testrunner_script_path.parent.parent.parent.resolve()
if "OGS_DATA_DIR" not in os.environ:
    os.environ["OGS_DATA_DIR"] = str(ogs_source_path / "Tests/Data")
Path(args.out).mkdir(parents=True, exist_ok=True)
build_dir = Path(args.out).parent.parent
success = True

for notebook_file in args.notebooks:
    notebook_file_path = Path(notebook_file)
    notebook_success = True
    is_jupytext = False
    if notebook_file_path.suffix in [".md", ".py"]:
        is_jupytext = True
    notebook_file_path_relative = notebook_file_path.absolute().relative_to(
        ogs_source_path
    )

    notebook_basename = notebook_file_path.parent.resolve() / notebook_file_path.stem
    _relpath = os.path.relpath(notebook_basename, start=os.environ["OGS_DATA_DIR"])
    notebook_output_path = (Path(args.out) / _relpath).resolve()
    notebook_output_path.mkdir(parents=True, exist_ok=True)
    os.environ["OGS_TESTRUNNER_OUT_DIR"] = str(notebook_output_path)
    os.environ["TQDM_DISABLE"] = "1"  # Disable progress bars
    os.environ["OGS_TESTRUNNER"] = "1"
    os.environ["PYVISTA_OFF_SCREEN"] = "true"
    os.environ.pop("OGS_TESTRUNNER_WEB_OUT_DIR", None)

    notebook_filename = notebook_file_path.name
    convert_notebook_file = notebook_output_path
    if not is_jupytext:
        convert_notebook_file = convert_notebook_file / Path(notebook_filename).stem
    convert_notebook_file = convert_notebook_file.with_suffix(".ipynb")

    nb = None
    executes_with_papermill = (
        not coverage_enabled() or notebook_file_path.suffix == ".md"
    )
    if "run-skip" not in str(notebook_file_path) and (
        args.hugo or executes_with_papermill
    ):
        if is_jupytext:
            nb = jupytext.read(notebook_file_path)
            convert_notebook_file = Path(
                str(convert_notebook_file).replace("notebook-", "")
            )
        else:
            with notebook_file_path.open(encoding="utf-8") as f:
                nb = nbformat.read(f, as_version=4)

        website_output_path = (
            get_website_output_path(nb, convert_notebook_file) if args.hugo else None
        )
        if args.hugo and website_output_path:
            os.environ["OGS_TESTRUNNER_WEB_OUT_DIR"] = str(website_output_path)

    if coverage_enabled() and notebook_file_path.suffix != ".md":
        os.environ["MPLBACKEND"] = "AGG"
        coverage_path = Path(build_dir / "coverage")
        coverage_path.mkdir(exist_ok=True)
        print(f"[Start]  {notebook_file}")
        run = subprocess.run(
            [
                "coverage",
                "run",
                f"--rcfile={ogs_source_path!s}/Tests/Data/pyproject.toml",
                f"--data-file={coverage_path.absolute().as_posix()}/.coverage",
                notebook_file,
            ],
            capture_output=True,
            cwd=notebook_file_path.parent,
            check=False,
        )

        if run.returncode == 0:
            print(f"[Passed]  {notebook_file}")
        else:
            print(f"[Failed] {notebook_file}.\n\n{run.stdout}\n\n{run.stderr}")
            sys.exit(1)

    elif "run-skip" not in str(notebook_file_path):
        if is_jupytext:
            jupytext.write(nb, convert_notebook_file)

        # Run the notebook
        print(f"[Start]  {notebook_filename}")
        start = timer()
        try:
            # Run with papermill instead of nbconvert for printing notebook
            # outputs on the command line
            nb = papermill.execute.execute_notebook(
                nb,
                None,
                kernel_name="python3",
                cwd=notebook_file_path.parent,
                log_output=True,
                progress_bar=False,
                stdout_file=sys.stdout,
                stderr_file=sys.stderr,
            )
        except DeadKernelError:
            out = None
            msg = f'Error executing the notebook "{notebook_filename}".\n\n'
            msg += f'See notebook "{convert_notebook_file}" for the traceback.'
            print(msg)
            notebook_success = False
            with convert_notebook_file.open(mode="w", encoding="utf-8") as f:
                nbformat.write(nb, f)
        except CellExecutionError:
            notebook_success = False
        end = timer()
        print(f"[End]  {notebook_filename}")

        # Write new notebook
        with convert_notebook_file.open(mode="w", encoding="utf-8") as f:
            if args.hugo:
                success = check_and_modify_frontmatter()
            nbformat.write(nb, f)

        status_string = ""
        if notebook_success:
            status_string += "[Passed] "
            if args.hugo:
                save_to_website(convert_notebook_file)
        else:
            status_string += "[Failed] "

            # Create and write HTML file
            html_exporter = HTMLExporter()
            html_exporter.template_name = "classic"
            body, resources = html_exporter.from_notebook_node(nb)

            html_file = convert_notebook_file.with_suffix(
                convert_notebook_file.suffix + ".html"
            )
            with html_file.open(mode="w", encoding="utf-8") as fh:
                fh.write(body)

        status_string += f"{notebook_filename} in "
        status_string += f"{timedelta(seconds=end-start).total_seconds()} seconds."
        if not notebook_success:
            status_string += f" --> {html_file} <--"
        print(status_string)


if not (success and notebook_success):
    sys.exit(1)
