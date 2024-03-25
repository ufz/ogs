+++
date = "{{ .Date }}"
title = "OpenGeoSys {{ .File.ContentBaseName }}"
tag = "{{ .File.ContentBaseName }}"
author = "Lars Bilke"
release_date = "{{ .Date }}"

[downloads]
binary = [
"Windows-10.0.22621-python-3.10.9-de-utils.zip",
"Windows-10.0.22621-python-3.10.9-utils.zip",
]
container = [
"serial.sif",
"openmpi-4.0.5.sif",
]
pip = true
note = """
**Note:** When using Python bindings make sure to have Python installed on your system:

- Windows: [Python 3.10.9](https://www.python.org/ftp/python/3.10.9/python-3.10.9-amd64.exe)
"""
+++
