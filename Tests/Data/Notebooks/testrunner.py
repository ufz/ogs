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


def save_to_website(exec_notebook_file):
    output_path_arg = ""
    output_path = ""
    notebook = nbformat.read(exec_notebook_file, as_version=4)
    first_cell = notebook.cells[0]
    if "Tests/Data" in str(exec_notebook_file):
        lines = first_cell.source.splitlines()
        toml_begin = lines.index("+++")
        toml_end = max(loc for loc, val in enumerate(lines) if val == "+++")
        toml_lines = lines[toml_begin + 1 : toml_end]
        parsed_frontmatter = toml.loads("\n".join(toml_lines))
        output_path = (
            Path(build_dir)
            / Path("web/content/docs/benchmarks")
            / Path(parsed_frontmatter["web_subsection"])
            / exec_notebook_file.stem.lower()
        )
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

    if "Tests/Data" not in str(exec_notebook_file):
        return

    Path(output_path).mkdir(parents=True, exist_ok=True)

    to_copy_path = Path(output_path) / "."

    if is_jupytext:
        shutil.copy(exec_notebook_file, to_copy_path)
        print(f"Copying ${exec_notebook_file} to {to_copy_path}")

    for subfolder in ["figures", "images"]:
        figures_path = (Path(notebook_file_path).parent / subfolder).resolve()
        symlink_figures_path = to_copy_path / subfolder
        if figures_path.exists() and not symlink_figures_path.exists():
            print(
                f"{subfolder} folder detected, copying {figures_path} to "
                f"{symlink_figures_path}"
            )
            shutil.copytree(figures_path, symlink_figures_path)


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
    first_cell.source = first_cell.source.replace("+++\n", "+++\nnotebook = true\n", 1)

    # Insert Jupyter header with notebook source and binderhub link in second cell
    repo = "https://gitlab.opengeosys.org/ogs/ogs"
    branch = "master"
    binder_tag = "6.5.6-0.7.0"
    if "CI_MERGE_REQUEST_SOURCE_PROJECT_URL" in os.environ:
        repo = os.environ["CI_MERGE_REQUEST_SOURCE_PROJECT_URL"]
        branch = os.environ["CI_MERGE_REQUEST_SOURCE_BRANCH_NAME"]
    binder_link = f"https://binder.opengeosys.org/v2/gh/bilke/binder-ogs-requirements/{binder_tag}?urlpath=git-pull%3Frepo={repo}%26urlpath=lab/tree/ogs/{notebook_file_path_relative}%26branch={branch}%26depth=1"
    text = """
<div class="note">
    <p style="margin-top: 0; margin-bottom: 0;">
        <img style="margin-top: 0; margin-bottom: 0; height: 2em;" class="inline-block mr-2 no-fancybox"
            src="https://upload.wikimedia.org/wikipedia/commons/3/38/Jupyter_logo.svg" alt="">
        This page is based on a Jupyter notebook."""
    if is_jupytext:
        download_file_name = Path(convert_notebook_file).name
        text += f"""
<a href="./{download_file_name}" download="{download_file_name}"><img class="no-fancybox" style="display: inline; margin-top: 0; margin-bottom: 0; margin-left: 1em;" src="https://img.shields.io/static/v1?label=Download:&message={download_file_name}&color=blue" /></a>"""
    text += f"""
<a href="{repo}/-/blob/{branch}/{notebook_file_path_relative}"><img src="https://img.shields.io/static/v1?label=Source:&message={notebook_filename}&color=brightgreen" class="no-fancybox"
        style="display: inline; margin-top: 0; margin-bottom: 0; margin-left: 1em;" /></a>
<a href="{binder_link}">
    <img class="no-fancybox" style="display: inline; margin-top: 0; margin-bottom: 0; margin-left: 1em;"
        src="https://img.shields.io/static/v1?label=&message=Launch notebook&color=5c5c5c&logo=data:image/png;base64,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" />
</a>"""
    text += """</p></div>\n\n"""

    second_cell = nb["cells"][1]
    if second_cell.cell_type == "markdown":
        second_cell.source = text + second_cell.source
    else:
        # Insert a new markdown cell
        nb["cells"].insert(1, nbformat.v4.new_markdown_cell(text))

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

    if "run-skip" not in str(notebook_file_path):
        notebook_basename = (
            notebook_file_path.parent.resolve() / notebook_file_path.stem
        )
        _relpath = os.path.relpath(notebook_basename, start=os.environ["OGS_DATA_DIR"])
        notebook_output_path = (Path(args.out) / _relpath).resolve()
        notebook_output_path.mkdir(parents=True, exist_ok=True)
        os.environ["OGS_TESTRUNNER_OUT_DIR"] = str(notebook_output_path)
        os.environ["TQDM_DISABLE"] = "1"  # Disable progress bars
        notebook_filename = notebook_file_path.name
        convert_notebook_file = notebook_output_path
        if not is_jupytext:
            convert_notebook_file = convert_notebook_file / Path(notebook_filename).stem
        convert_notebook_file = convert_notebook_file.with_suffix(".ipynb")

        if is_jupytext:
            nb = jupytext.read(notebook_file_path)
            convert_notebook_file = Path(
                str(convert_notebook_file).replace("notebook-", "")
            )
            jupytext.write(nb, convert_notebook_file)
        else:
            with notebook_file_path.open(encoding="utf-8") as f:
                nb = nbformat.read(f, as_version=4)

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
        (body, resources) = html_exporter.from_notebook_node(nb)

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
