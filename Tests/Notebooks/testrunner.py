import nbformat
from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError
from nbconvert import HTMLExporter
import argparse
import os
import sys
from timeit import default_timer as timer
from datetime import timedelta


# Script arguments
parser = argparse.ArgumentParser(description="Jupyter notebook testrunner.")
parser.add_argument("notebooks", metavar="N", nargs="+", help="Notebooks to test.")
parser.add_argument("--out", default="./", help="Output directory.")
parser.add_argument(
    "--hugo", action="store_true", help="Convert successful notebooks to web site."
)
args = parser.parse_args()

# Path setup
testrunner_script_path = os.path.dirname(os.path.abspath(__file__))
ogs_source_path = os.path.abspath(testrunner_script_path + "/../..")
if "OGS_DATA_DIR" not in os.environ:
    os.environ["OGS_DATA_DIR"] = ogs_source_path + "/Tests/Data"
os.environ["OGS_TESTRUNNER_OUT_DIR"] = args.out
os.makedirs(args.out, exist_ok=True)
success = True

for notebook_file_path in args.notebooks:
    notebook_success = True

    with open(notebook_file_path) as f:
        nb = nbformat.read(f, as_version=4)
    ep = ExecutePreprocessor(timeout=600, kernel_name="python3")

    # 1. Run the notebook
    start = timer()
    try:
        ep.preprocess(nb, {"metadata": {"path": "."}})
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

    # write new notebook
    # TODO preserve relative file structure in output
    notebook_filename = os.path.basename(notebook_file_path)
    exec_notebook_file = args.out + "/" + notebook_filename
    with open(exec_notebook_file, "w", encoding="utf-8") as f:
        nbformat.write(nb, f)

    status_string = ""
    if notebook_success:
        status_string += "[Passed] "
        if args.hugo:
            # Save to website
            from nb2hugo.writer import HugoWriter

            writer = HugoWriter()
            web_path = ogs_source_path + "/web"
            writer.convert(
                exec_notebook_file,
                web_path,
                "docs/benchmarks/notebooks",
                testrunner_script_path + "/nbconvert_templates/collapsed.md.j2",
            )
            # TODO do not write if it has no frontmatter OR return frontmatter with draft: true
    else:
        status_string += "[Failed] "
        # Save to html
        html_file = exec_notebook_file + ".html"
        with open(html_file, "w", encoding="utf-8") as fh:
            fh.write(body)

    status_string += f"{notebook_filename} in {timedelta(seconds=end-start).total_seconds()} seconds."
    if not notebook_success:
        status_string += f" --> {html_file} <--"
    print(status_string)

if not success:
    sys.exit(1)
