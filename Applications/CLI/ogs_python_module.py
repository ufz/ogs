#!/usr/bin/env python3

# The OpenGeoSys python library will be created in the lib folder of the build
# directory
## export PYTHONPATH=$PYTHONPATH:your-build_directory-here/lib

import sys
import tempfile
import OpenGeoSys

arguments = ["", sys.argv[1], "-o " + tempfile.mkdtemp()]

print("Python OpenGeoSys.init ...")
OpenGeoSys.initialize(arguments)
print("Python OpenGeoSys.executeSimulation ...")
OpenGeoSys.executeSimulation()
print("Python OpenGeoSys.finalize() ...")
OpenGeoSys.finalize()
print("Python world.")
