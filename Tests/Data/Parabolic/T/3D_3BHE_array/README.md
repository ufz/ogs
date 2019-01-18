# How to run

You need to build OGS with `OGS_USE_PYTHON=ON` and have the Python modules
listed in `requirements.txt` installed, e.g.:

- In this directory:
  - `virtualenv .venv` Make sure you use the same Python interpreter which was
    embedded into OGS (look for `Found PythonInterp: ..` in the CMake output)
  - `.venv/bin/activate`
  - `pip install -r requirements.txt`
  - `PYTHONPATH=.venv/lib/python3.7/site-packages ~/code/ogs6/build-n/bin/ogs -o ~/tmp 3bhes.prj`
