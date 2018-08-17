+++
project = "Parabolic/ComponentTransport/elder/elder-python.prj"
author = "Christoph Lehmann"
date = "2018-08-16T09:18:00+02:00"
title = "Saturated Variable-Density Flow and Mass Transport (Elder) with Python BC"
weight = 3

[menu]
  [menu.benchmarks]
    parent = "python-bc"

+++

{{< data-link >}}

## Motivation of this test case

The aim of this test is:

* to show that it is possible to prescribe BCs only on parts of a given geometry
* to assert that Python BCs work with processes involving more than one physical
  field.

## Details

This test is a copy of [this test case]({{< ref "../../hydro-component/elder.md" >}}).
Please check the original test case for any details.
