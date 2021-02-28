+++
title = "Add element quality to mesh"
author = "Thomas Fischer"

[menu]
  [menu.tools]
    parent = "model-preparation"
+++

The command line tool `AddElementQuality` adds a data array to the given mesh
that holds for each element a value for the quality. The quality value depends
on the chosen quality criterion. The available quality criteria are shown using
the help argument on the command line. The resulting data arrays can be used to
draw histograms or can be visually analyzed in the OGS Data Explorer or Paraview.
