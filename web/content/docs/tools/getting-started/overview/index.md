+++
date = "2018-02-27T11:00:13+01:00"
title = "Introduction"
author = "Lars Bilke"
weight = 1
aliases = ["/docs/tools/"]

[menu.docs]
name = "Tools & Workflows"
identifier = "tools"
weight = 4
post = "Helpful tools for pre- and postprocessing as well as complete model setup workflows."
[menu.docs.params]
category = "User"
+++

Here is an overview of the currently available command line tools, that might help you to construct your OpenGeoSys model. GUI-based tools are available too:

- OpenGeoSys Data Explorer -- [Manual](https://gitlab.opengeosys.org/ogs/documentation/data_explorer_manual/-/jobs/artifacts/master/raw/ogsde-man.pdf?job=build) / [Download](/releases)
- [GINA by BGR](https://teambeam.bgr.de/my/drive/folder/68)

----

You'll find examples and applications of the tools. Choose the appropriate tool to the left. Ready to use binaries can be downloaded on the [Releases](/releases)-page. The second possibility to obtain the tools is to check out OGS sources and compile the tools.

Feel free to get in touch with us if you have any issues with any tool.

The tools are ordered similarly to a workflow to setup might look.

At the beginning of the modeling often external data have to be converted into OGS data formats / OGS data structures. Tools helping the modeler are listed in the data conversion tools section.

This data can be used to create meshes approximating geometrically the domain of interest. Despite some simple mesh creator tools OGS offers some functionality to analyze meshes. Furthermore there are some tools for cleanup the mesh, adding a layer at the top or moving the mesh.

The next important step creating a simulation model is the parametrization of the model, i.e., material parameters, boundary conditions and source terms have to be assigned to the model.
