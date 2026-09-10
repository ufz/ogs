+++
title = 'OpenWorkFlow'
abbreviation = 'OWF'
headline = 'Synthesis platform for the German site selection'

links = [
  ['OpenWorkFlow project page', 'https://www.openworkflow.de/'],
]
+++

According to the German Site Selection Act (StandAG), the objective is to identify the site that provides the highest possible level of safety for a repository for high-level radioactive waste in Germany. To this end, safety assessments are conducted for a period of one million years, among other measures. Many of the analyses required for these assessments, in turn, rely on numerical simulations of the physical, chemical, and biological processes occurring within a repository and its surroundings.

An important project objective is the development of prototype models for potential repositories/site regions. These consist of linked numerical models whose key outputs are safety-relevant quantities, such as material/mass releases, fluid pressure, dilatancy, and temperature. The prototype models differ for the individual host rocks—clay rock, rock salt, and crystalline rock. They serve as blueprints and can be readily adapted to specific site regions. To achieve this project objective, OpenWorkFlow is developing a numerical toolbox for safety assessments that covers tasks such as data integration, model creation, model execution, post-processing, and archiving.

Both the development process and its outcome (the synthesis platform) encompass a broad range of scientific and technical aspects and are subject to a number of quality criteria, which are described in more detail in the following sections.

More information can be found on the [OpenWorkFlow web site](https://www.openworkflow.de) (in German only).

{{< figure src="owf-overview.jpg" caption="Schematic representation of key aspects of the developed simulation workflows, shown clockwise. *a*: Integration of (geological) data into computational models. *b*: Simulation of FEPs (here: glaciation). *c*: Coupled thermal, hydraulic, mechanical, and chemical processes in the far field. *d*: Relevant geometries in the near field—fuel assemblies, canisters, host rock, etc. *e*: Analysis of geotechnical support systems—grouted anchors. *f*: Near-field integrity analysis—gas pressure increase due to canister corrosion for different parameter sets. *g*: Tracing the execution graph of a workflow using the provenance database." >}}

## References

- {{< bib "Heinze2025" >}}
- {{< bib "buchwald2025" >}}
- {{< bib "buchwald2024" >}}
- {{< bib "Lehmann2024" >}}
- {{< bib "Selzer2024" >}}
