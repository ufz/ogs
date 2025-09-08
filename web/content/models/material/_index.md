+++
title = 'Materials'
abbreviation= 'MAT'
links = [
  ['Hoek-Brown yield criterion', 'https://www.opengeosys.org/docs/benchmarks/small-deformations/hoekbrownyieldcriterion/'],
  ['TODO: Mohr-Coulomb', '#'],
  ['TODO: Cam-Clay', '#']
]
+++

Clay - Salt - Crystalline. <span class="text-xs">Image by BGE</span>

<div class="not-prose">
{{< figure
  src="mat.png"
  class="float-right"
>}}
</div>

Constitutive models for porous and fractured media describe the material behavior of their solid and fluid phases. Material models close the balance equation to thermodynamically consistent models, are independent of the corresponding conservation principles, and can therefore be implemented separately. This is done, for example, in [MFront](https://www.sciencedirect.com/science/article/pii/S0898122115003132?via%3Dihub) and the MPL library from OGS. Examples of tests of the constitutive models are:
