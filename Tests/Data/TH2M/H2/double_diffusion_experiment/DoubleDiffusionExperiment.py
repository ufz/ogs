# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.1
#   kernelspec:
#     display_name: python_3.11
#     language: python
#     name: python3
# ---

# %% [raw]
# +++
# title = "Double diffusion experiment in Boom Clay"
# date = "2026-08-27"
# author = "Norbert Grunwald, Michael Putz"
# image = "figures/double_diffusion_setup.png"
# web_subsection = "th2m"
# projects = ["eurad"]
# models = ["lab"]
# +++


# %% [markdown]
"""
<div class="note">

### Note: Placeholder benchmark

This notebook documents the planned OpenGeoSys representation of the double
through-diffusion experiment. The original project files and reference results
will be added when supplied by the authors. It intentionally does not run an
OGS simulation yet.

</div>

Two gas species diffuse in opposite directions through a fully water-saturated
Boom Clay sample. Helium is initially supplied from the left reservoir and
methane from the right reservoir. The experiment was performed at SCK CEN and
reported by Pitz et al. (2024).

![Schematic of the double diffusion experiment](./figures/double_diffusion_setup.svg "Double diffusion experiment: helium and methane diffuse through a saturated Boom Clay core in opposite directions.")

## Project context

The benchmark was developed within the European Joint Programme on Radioactive
Waste Management (EURAD), Work Package GAS (EURAD WP-GAS). The project was
funded by the Horizon 2020 Euratom programme under grant agreement No. 847593
(2019--2024).

## Experiment at a glance

| Item | Description |
| --- | --- |
| Host material | Fully water-saturated Boom Clay |
| Specimen | Cylinder, length 30 mm, radius 40 mm |
| Reservoirs | Two 1 L vessels with water and a gas phase |
| Initial gases | Helium (left) and methane (right), initially 1 MPa total pressure |
| Transport | Dissolution in pore water and molecular diffusion; no imposed pressure-driven advection |
| Observations | Time-dependent reservoir pressure and gas composition, measured by gas chromatography over 72 days |

Under the experiment's conditions, each gas moves down its own concentration
gradient: helium towards the methane vessel and methane towards the helium
vessel. Gas--water partitioning is assumed to be at instantaneous equilibrium
and is described by Henry's law. The downstream gas concentration is obtained
from the accumulated component mass leaving the clay sample, together with the
water and gas volumes of that reservoir.

## Reference

{{< bib "Pitz2024" >}}

"""
