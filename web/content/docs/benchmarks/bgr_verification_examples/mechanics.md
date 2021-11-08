+++
date = "2019-04-24T15:15:45+01:00"
title = "Small deformations: Verification examples by Vogel,Maßmann"
weight = 50
author = "Johannes Herfurth, Jan Thiedau"

[menu]
  [menu.benchmarks]
    parent = "small-deformations"

+++

These benchmark examples test the implementation of
small deformations process with analytical solutions
presented by Vogel/Massmann.

A detailed description can be found in the ogs Benchmark books.
The following table links the ogs problem descriptions with its corresponding
chapters in the benchmark books.

| Book/Chapter | Benchmark name |
|:--- | :--- |
|*Kolditz et al. 2015*||
|2.4.1 | m1_1Dload|
|2.4.3 | m1_3Dgravity|
| *Kolditz et al. 2016*||
|4.1.1   |  m1_1Dlozenge|
|4.1.2   |  m1_2Dload|
|4.1.3   |  m1_3Dload|
| *Kolditz et al. 2018*||
|3.7  |  m1_3Dsquare|
|3.9  |  m1_3Dbottom|
|3.10 |  m1_3Dtopload|
|3.15  |           {{< data-link "m3_3Dshearz" "Mechanics/Linear/Orthotropy/m3_3Dshearz.prj" >}}|
|3.15 (rotated) |  {{< data-link "m3_3Dshearz_rot" "Mechanics/Linear/Orthotropy/m3_3Dshearz_rot.prj" >}}|
|3.16  |           {{< data-link "m3_3Dtopload" "Mechanics/Linear/Orthotropy/m3_3Dtopload.prj" >}}|
|3.16 (rotated) |  {{< data-link "m3_3Dtoploadlc" "Mechanics/Linear/Orthotropy/m3_3Dtoploadlc.prj" >}}|

Note that OGS uses plane strain assumptions in two dimensions (2D) as it is common in geomechanics.

## References

{{< bib "kolditz:2015" >}}

{{< bib "kolditz:2016" >}}

{{< bib "kolditz:2018" >}}
