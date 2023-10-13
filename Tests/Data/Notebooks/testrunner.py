import nbformat
from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError
from nbclient.exceptions import DeadKernelError
from nbconvert import HTMLExporter
import argparse
import os
import shutil
import sys
from timeit import default_timer as timer
from datetime import timedelta
import toml
from pathlib import Path
import jupytext
import subprocess


def save_to_website(exec_notebook_file, web_path):
    output_path_arg = ""
    notebook = nbformat.read(exec_notebook_file, as_version=4)
    first_cell = notebook.cells[0]
    if "Tests/Data" in exec_notebook_file:
        lines = first_cell.source.splitlines()
        toml_begin = lines.index("+++")
        toml_end = max(loc for loc, val in enumerate(lines) if val == "+++")
        toml_lines = lines[toml_begin + 1 : toml_end]
        parsed_frontmatter = toml.loads("\n".join(toml_lines))
        output_path = (
            Path(build_dir)
            / Path("web/content/docs/benchmarks")
            / Path(parsed_frontmatter["web_subsection"])
        )
        output_path_arg = (
            f"--output-dir={Path(output_path) / Path(exec_notebook_file).stem}"
        )

    template = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "nbconvert_templates/collapsed.md.j2",
    )
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

    if not "Tests/Data" in exec_notebook_file:
        return

    for subfolder in ["figures", "images"]:
        figures_path = os.path.abspath(
            os.path.join(os.path.dirname(notebook_file_path), subfolder)
        )
        symlink_figures_path = os.path.join(
            web_path,
            "content",
            output_path,
            os.path.splitext(os.path.basename(exec_notebook_file))[0].lower(),
            subfolder,
        )
        if os.path.exists(figures_path) and not os.path.exists(symlink_figures_path):
            print(
                f"{subfolder} folder detected, copying {figures_path} to "
                f"{symlink_figures_path}"
            )
            shutil.copytree(figures_path, symlink_figures_path)


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
testrunner_script_path = os.path.dirname(os.path.abspath(__file__))
ogs_source_path = os.path.abspath(os.path.join(testrunner_script_path, "../../.."))
if "OGS_DATA_DIR" not in os.environ:
    os.environ["OGS_DATA_DIR"] = os.path.join(ogs_source_path, "Tests/Data")
os.makedirs(args.out, exist_ok=True)
build_dir = Path(args.out).parent.parent
success = True

for notebook_file_path in args.notebooks:
    notebook_success = True
    is_jupytext = False
    if Path(notebook_file_path).suffix in [".md", ".py"]:
        is_jupytext = True
    notebook_file_path_relative = (
        Path(notebook_file_path).absolute().relative_to(ogs_source_path)
    )

    if "run-skip" not in notebook_file_path:
        notebook_basename = os.path.splitext(notebook_file_path)[0]
        notebook_output_path = os.path.abspath(
            os.path.join(
                args.out,
                os.path.relpath(notebook_basename, start=os.environ["OGS_DATA_DIR"]),
            )
        )
        os.makedirs(notebook_output_path, exist_ok=True)
        os.environ["OGS_TESTRUNNER_OUT_DIR"] = notebook_output_path
        notebook_filename = os.path.basename(notebook_file_path)
        convert_notebook_file = notebook_output_path
        if not is_jupytext:
            convert_notebook_file = os.path.join(
                convert_notebook_file, Path(notebook_filename).stem
            )
        convert_notebook_file += ".ipynb"

        if is_jupytext:
            nb = jupytext.read(notebook_file_path)
            convert_notebook_file = convert_notebook_file.replace("notebook-", "")
            jupytext.write(nb, convert_notebook_file)
        else:
            with open(notebook_file_path, mode="r", encoding="utf-8") as f:
                nb = nbformat.read(f, as_version=4)
        ep = ExecutePreprocessor(kernel_name="python3")

        # Run the notebook
        print(f"[Start]  {notebook_filename}")
        start = timer()
        try:
            ep.preprocess(
                nb, {"metadata": {"path": os.path.dirname(notebook_file_path)}}
            )
        except DeadKernelError:
            out = None
            msg = 'Error executing the notebook "%s".\n\n' % notebook_filename
            msg += 'See notebook "%s" for the traceback.' % convert_notebook_file
            print(msg)
            notebook_success = False
            success = False
            with open(convert_notebook_file, mode="w", encoding="utf-8") as f:
                nbformat.write(nb, f)
        except CellExecutionError:
            notebook_success = False
            success = False
            pass
        end = timer()

        # Write new notebook
        with open(convert_notebook_file, "w", encoding="utf-8") as f:
            repo = "https://gitlab.opengeosys.org/ogs/ogs"
            branch = "master"
            if "CI_MERGE_REQUEST_SOURCE_PROJECT_URL" in os.environ:
                repo = os.environ["CI_MERGE_REQUEST_SOURCE_PROJECT_URL"]
                branch = os.environ["CI_MERGE_REQUEST_SOURCE_BRANCH_NAME"]

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

            # Check second cell is markdown
            second_cell = nb["cells"][1]
            if second_cell.cell_type != "markdown":
                print(
                    f"Error: {notebook_filename} first cell after the frontmatter needs to be a markdown cell! Move the first Python cell below."
                )
                success = False

            # Modify metadata
            first_cell = nb["cells"][0]
            if first_cell.source.startswith("---"):
                print(
                    f"Error: {notebook_filename} frontmatter is not in TOML format! Use +++ delimitiers!"
                )
                success = False
            first_cell.source = first_cell.source.replace(
                "+++\n", "+++\nnotebook = true\n", 1
            )

            # Insert Jupyter header with notebook source and binderhub link
            binder_link = f"https://mybinder.org/v2/gh/bilke/binder-ogs-requirements/master?urlpath=git-pull%3Frepo={repo}%26urlpath=lab/tree/ogs/{notebook_file_path_relative}%26branch={branch}"
            text = f"""
<div class="note">
    <p style="margin-top: 0; margin-bottom: 0;">
        <img style="margin-top: 0; margin-bottom: 0; height: 2em;" class="inline-block mr-2 no-fancybox"
            src="https://upload.wikimedia.org/wikipedia/commons/3/38/Jupyter_logo.svg" alt="">
        This page is based on a Jupyter notebook."""
            if is_jupytext:
                download_file_name = (
                    Path(convert_notebook_file)
                    .rename(Path(convert_notebook_file).with_suffix(".ipynb"))
                    .name
                )
                text += f"""
<a href="./{download_file_name}" download="{download_file_name}"><img class="no-fancybox" style="display: inline; margin-top: 0; margin-bottom: 0; margin-left: 1em;" src="https://img.shields.io/static/v1?label=Download:&message={download_file_name}&color=blue" /></a>"""
            text += f"""
<a href="{repo}/-/blob/{branch}/{notebook_file_path_relative}"><img src="https://img.shields.io/static/v1?label=Source:&message={notebook_filename}&color=brightgreen" class="no-fancybox"
        style="display: inline; margin-top: 0; margin-bottom: 0; margin-left: 1em;" /></a>
<a href="{binder_link}">
    <img class="no-fancybox" style="display: inline; margin-top: 0; margin-bottom: 0; margin-left: 1em;"
        src="https://img.shields.io/static/v1?label=&message=Launch notebook&color=5c5c5c&logo=data:image/png;base64,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" />
</a>"""
            text += f"""</p></div>\n\n"""
            second_cell.source = text + second_cell.source

            nbformat.write(nb, f)

    status_string = ""
    if notebook_success:
        status_string += "[Passed] "
        if args.hugo:
            save_to_website(
                convert_notebook_file, os.path.join(ogs_source_path, args.hugo_out)
            )
    else:
        status_string += "[Failed] "
        # Save to html
        html_file = convert_notebook_file + ".html"
        with open(html_file, "w", encoding="utf-8") as fh:
            fh.write(body)

    status_string += f"{notebook_filename} in "
    status_string += f"{timedelta(seconds=end-start).total_seconds()} seconds."
    if not notebook_success:
        status_string += f" --> {html_file} <--"
    print(status_string)

if not success:
    sys.exit(1)
