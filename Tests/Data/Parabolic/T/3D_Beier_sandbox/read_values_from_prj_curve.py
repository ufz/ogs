# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import inspect
import os
import xml.etree.ElementTree as ET

import numpy as np

file = inspect.getfile(inspect.currentframe())
prj_path = os.path.split(os.path.realpath(file))[0]

# read prj-data
tree = ET.parse(prj_path + "/beier_sandbox.prj")
root = tree.getroot()

curve_elements = root.findall(".//curve")
for curve in curve_elements:
    curve_name = curve.find("name").text

    # read coordinates
    coords_element = curve.find("coords")
    coords_text = coords_element.text  # replace("\n", "    ")
    coords_array = np.fromstring(coords_text, dtype=float, sep="    ")
    # Store array in binary format
    coords_array.astype("<f8").tofile(
        prj_path + f"/{curve_name}_coords.bin"
    )  # < - little Endian | f - float | 8 - Bytes equal to 64bit(double)

    # read values
    values_element = curve.find("values")
    values_text = values_element.text  # replace("\n", "    ")
    values_array = np.fromstring(values_text, dtype=float, sep="    ")
    values_array.astype("<f8").tofile(
        prj_path + f"/{curve_name}_values.bin"
    )  # < - little Endian | f - float | 8 - Bytes equal to 64bit(double)
