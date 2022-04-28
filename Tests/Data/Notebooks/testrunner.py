import nbformat
from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError
from nbconvert import HTMLExporter
import argparse
import os
import sys
from timeit import default_timer as timer
from datetime import timedelta


def save_to_website(exec_notebook_file, web_path):
    from nb2hugo.writer import HugoWriter

    # TODO do not write if it has no frontmatter OR return frontmatter with draft: true
    writer = HugoWriter()
    writer.convert(
        exec_notebook_file,
        web_path,
        "docs/benchmarks/notebooks",
        os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "nbconvert_templates/collapsed.md.j2",
        ),
    )


# Script arguments
parser = argparse.ArgumentParser(description="Jupyter notebook testrunner.")
parser.add_argument("notebooks", metavar="N", nargs="+", help="Notebooks to test.")
parser.add_argument("--out", default="./", help="Output directory.")
parser.add_argument(
    "--hugo", action="store_true", help="Convert successful notebooks to web site."
)
parser.add_argument(
    "--timeout", default="600", type=int, help="Cell execution timeout."
)
args = parser.parse_args()

# Path setup
testrunner_script_path = os.path.dirname(os.path.abspath(__file__))
ogs_source_path = os.path.abspath(os.path.join(testrunner_script_path, "../../.."))
if "OGS_DATA_DIR" not in os.environ:
    os.environ["OGS_DATA_DIR"] = os.path.join(ogs_source_path, "Tests/Data")
os.makedirs(args.out, exist_ok=True)
success = True

for notebook_file_path in args.notebooks:
    notebook_success = True
    convert_notebook_file = notebook_file_path

    if "run-skip" not in notebook_file_path:
        notebook_basename = os.path.splitext(notebook_file_path)[0]
        notebook_output_path = os.path.join(
            args.out,
            os.path.relpath(notebook_basename, start=os.environ["OGS_DATA_DIR"]),
        )
        os.makedirs(notebook_output_path, exist_ok=True)
        os.environ["OGS_TESTRUNNER_OUT_DIR"] = notebook_output_path

        with open(notebook_file_path) as f:
            nb = nbformat.read(f, as_version=4)
        ep = ExecutePreprocessor(timeout=args.timeout, kernel_name="python3")

        # 1. Run the notebook
        start = timer()
        try:
            ep.preprocess(
                nb, {"metadata": {"path": os.path.dirname(notebook_file_path)}}
            )
        except CellExecutionError:
            notebook_success = False
            success = False
            pass
        end = timer()

        # 2. Instantiate the exporter. We use the `classic` template for now; we'll get into more details
        # later about how to customize the exporter further.
        html_exporter = HTMLExporter()
        html_exporter.template_name = "classic"

        # 3. Process the notebook we loaded earlier
        (body, resources) = html_exporter.from_notebook_node(nb)

        # 4. Write new notebook
        notebook_filename = os.path.basename(notebook_file_path)
        convert_notebook_file = os.path.join(notebook_output_path, notebook_filename)
        with open(convert_notebook_file, "w", encoding="utf-8") as f:
            nbformat.write(nb, f)

    status_string = ""
    if notebook_success:
        status_string += "[Passed] "
        if args.hugo:
            save_to_website(convert_notebook_file, os.path.join(ogs_source_path, "web"))
    else:
        status_string += "[Failed] "
        # Save to html
        html_file = convert_notebook_file + ".html"
        with open(html_file, "w", encoding="utf-8") as fh:
            fh.write(body)

    status_string += f"{notebook_filename} in {timedelta(seconds=end-start).total_seconds()} seconds."
    if not notebook_success:
        status_string += f" --> {html_file} <--"
    print(status_string)

if not success:
    sys.exit(1)
