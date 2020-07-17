+++
date = "2019-12-03T00:00:00+01:00"
title = "Mesh2Shape"
author = "Karsten Rink"

[menu]
  [menu.tools]
    parent = "Data Import/Export"
+++

## Introduction

Converts a 2D surface mesh into a shapfile such that each element is represented by a polygon. Cell attributes are transferred onto shape polygons while point attributes are ignored.

## Usage

```bash
   Mesh2Shape  -i <input_file.vtu> -o <output_file.shp>


Where:

   -i <input_file.vtu>,  --input-file <input_file.vtu>
     (required)  OGS mesh file (*.vtu, *.msh)

   -o <output_file.shp>,  --output-file <output_file.shp>
     (required)  Esri Shapefile (*.shp)
```

## Simple example

**Input data:**

![2D surface mesh with scalar data assigned to cells, here displayed via the OGS Data Explorer. In this particular case, the simulation result of groundwater flow simulation (originally assigned to mesh nodes) has been converted onto cells via VTK's PointToCell-Filter.](./Mesh2Shape-input.png){.m-auto}

**Command:**

```bash
Mesh2Shape -i Mueglitz2D_Point2Cell.vtu -o Mueglitz2D_Point2Cell.shp
```

![Exported shapefile displayed in a geographic information system (here, QGIS).](./Mesh2Shape-output1.png){.m-auto}

![The result of an OGS-simulation showing the groundwater head of the Müglitz-catchment imported into QGIS and combined with other data from an existing GIS-project of this region.](./Mesh2Shape-output2.png){.m-auto}

## Application

The utility allows to export meshes, and in particular simulation results, into existing GIS-projects and use the data as input for subsequent workflows.
