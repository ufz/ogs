+++
date = "2020-06-22T11:00:13+01:00"
title = "AddLayer"
author = "Thomas Fischer"
aliases = [ "/docs/tools/meshing/addtoplayer" ]
+++

## Introduction

The tool `AddLayer` adds one layer of elements with a specified thickness
`thickness` on either on top or bottom of an existing mesh `input-mesh` and
returns the newly generated mesh `output-mesh` that has an new layer.

One might want to take care that the material groups are reduced, e.g. material
groups should not be [0,2,5], but [0,1,2]. The new layer will have the material
group ID of the highest material group +1 if not otherwise specified. With
the switch `--copy-material-ids` the material group IDs of the extruded layer will be
kept. With the switch `--set-material-id` the material group id of the new layer can be set
to a different value than the default value of 'highest material group number plus one'.

The tool requires the [OGS-6 node ordering](http://doxygen.opengeosys.org/index.html) in the elements. A different node ordering may lead to unexpected results. In case one might try to change the ordering using the tool `NodeReordering`.

## Usage

```bash
AddLayer -i <input-mesh> -o <output-mesh> [-t <thickness>]
```

## Simple example

![One material](SimpleQuadExample_1.png#one-third "A simple cube mesh with one material group (red).")

![Added layer](SimpleQuadExampleWithNewTopLayer_1.png#one-third "The updated mesh where an additional layer (blue) was added on top of the domain with a second material group.")

Usage for example:

```bash
AddLayer -i quad.vtu -o quad_with_new_top_layer.vtu -t 1
```

## Application

The tool was used to add a "soil" layer to the hydro-geological model of the Unstrut catchment within the INFLUINS project:

{{< bib "fischer:2015" >}}

The tool is used in the OGS-GIScape workflow.
