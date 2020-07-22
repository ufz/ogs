+++
date = "2019-12-03T00:00:00+01:00"
title = "createIntermediateRasters"
author = "Karsten Rink"

[menu]
  [menu.tools]
    parent = "Preprocessing"
+++

## Introduction

This utility uses two input DEMs (i.e. raster files of digital elevation models) located at the exact same spatial position but at different elevation and calculates a specified number of raster DEMs located at equidistant distance between them (i.e. for n=1, one new raster located precisely in the middle will be created).

## Usage

```bash
   createIntermediateRasters --file1 <file1.asc> --file2 <file2.asc>
                             -o <output.asc> [-n <int>]

Where:

   --file1 <file1.asc>
     (required)  First DEM-raster file

   --file2 <file2.asc>
     (required)  Second DEM-raster file

   -o <output.asc>,  --output-file <output.asc>
     (required)  Raster output file (*.asc)

   -n <int>,  --number <int>
     number of rasters to be calculated
```

The parameter ```n``` determines how many layers are created between the two input layers. If the parameter is not specified the default value is set as ```n=1```.

## Simple example

**Input data:**

![Two input rasters as well as their 3D surface representation. Darker pixels represent values at a lower elevation while brighter pixels represent higher elevaton. In the 3D visualisation, the left-most raster is represented by the green surface and the right-most raster by the yellow surface.](./createIntermediateRasters-input.png)

**Command:**

```bash
createIntermediateRasters --file1 raster1.asc --file2 raster2.asc -o output.asc -n 1
```

![A new raster is created in the exact center between the two input rasters. In the 3D representation, the new layer is shown in red.](./createIntermediateRasters-output1.png){ width=66% }

**Command:**

```bash
createIntermediateRasters --file1 raster1.asc --file2 raster2.asc -o output.asc -n 2
```

![For ```n>1``` multiple rasters are created at equidistant distances between the two input rasters. For ```n=2```, two new rasters are generated, represented here in red and blue.](./createIntermediateRasters-output2.png)

## Application

This utility allows to generate additional input data when creating a bulk mesh from geometry and raster layers. The existing layer structure given via raster files can be refined to avoid deformed elements if the resolution in x-/y-direction is significantly more detailed than in z-direction.
