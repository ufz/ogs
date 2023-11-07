+++
date = "2018-03-07T15:59:57+01:00"
title = "Remove Mesh Elements"
author = "Thomas Fischer"
+++

## General

The tool `removeMeshElements` removes those elements from a given input mesh that fulfills a user specified criterion. The resulting mesh will be written to the specified output file. The user can choose between 4 different removal criteria:

 1. Remove elements based on assigned properties, for instance material ids. Elements with properties outside of the given range are removed.
 2. Remove elements by element type, for instance remove line elements.
 3. Remove elements that have zero volume.
 4. Remove elements based on an axis aligned bounding box. Elements where at least one point is located outside the specified bounding box (or - if the "invert"-flag is set - inside the bounding box) are removed.

One possible application is to cut out a smaller mesh out of a bigger one either by specifying a bounding box or by marking the inner/outer region with a unique MaterialID using the tool [SetPropertiesInPolygonalRegion]({{< ref "set-properties-in-polygonal-region" >}}).

Another application is to cut out patches of a (top) surface (tool [ExtractSurface]({{< ref "extract-surface" >}})) for assigning boundary conditions.

## Usage

```bash
removeMeshElements -i <input-mesh> -o <output-mesh>
 [-n <property_name>] [--int-property-value <number value>] ...
 [--element-type <element type>] ...
 [--zero-volume] ...
 [--x-min <value>] [--x-max <value>] [--y-min <value>] [--y-max <value>] [--z-min <value>] [--z-max <value>] [--invert]
```

Each particular line with optional arguments refers to one of the different removal criteria mentioned in the general section.
The corresponding element types differ from VTK cell types and can be found in `MeshLib/MeshEnums.cpp`.

## Examples

![Input](ExampleRemoveElements-Input.png "The left figure above is the result of the repeated application of [SetPropertiesInPolygonalRegion]({{< ref "set-properties-in-polygonal-region" >}}). It contains material ids 0 (red), 1 (yellow), 2 (turquoise) and 3 (blue).")

![The result of the following command line input is depicted.](ExampleRemoveElements-Output.png "The result of the following command line input is depicted.")

```bash
removeMeshElements -i TestCube-ResetPropertiesInPolygonalRegion.vtu -o TestCube-removeMeshElements.vtu -n MaterialIDs --int-property-value 1 --int-property-value 2 --int-property-value 3
```

## Applications

The tool was used to cut the Unstrut catchment out of the Thuringian syncline model and to remove some geological units not necessary for the modeling within the INFLUINS project, see reference [GO2OGS].

<div class='note'>

### Example Files

[TestCube-removeMeshElements.vtu](TestCube-removeMeshElements.vtu)  
[TestCube-ResetPropertiesInPolygonalRegion](TestCube-ResetPropertiesInPolygonalRegion.vtu)  
</div>
