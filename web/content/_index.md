---
title: OpenGeoSys
authors:
  - Lars Bilke

date: 2017-01-13T14:24:23+01:00

hero:
  eyebrow: Open-source multi-physics · OGS-6
  headline: Simulating coupled processes in porous & fractured media
  textline: OpenGeoSys (OGS) is a scientific [open-source](https://gitlab.opengeosys.org/ogs/ogs) finite-element framework for the numerical simulation of thermo-hydro-mechanical-chemical (THMC) processes in the subsurface.
  images:
    - permalink: features/vis/layeredview2.png
      alt: Layered geological model
      caption: layered geological model · ParaView / VTK
    - permalink: features/vis/DEMSelke3D.png
      alt: Discrete element groundwater model
      caption: discrete element model · groundwater domain
    - permalink: features/vis/contours2-bw.png
      alt: Simulation contour visualisation
      caption: contour visualisation of simulation results
  cta:
    - text: Download OGS-6
      url: /releases
    - text: Read the docs
      url: /docs

stats:
  - value: BSD-3-Clause
    label: Free & open source
    link: https://gitlab.opengeosys.org/ogs/ogs

application_areas_intro: The same THMC framework underpins research questions across the geosciences — three of the most active fields are shown below.

application_areas:
  - eyebrow: 01 / Subsurface water
    headline: Hydrology
    textline: From regional groundwater flow to contaminant and coastal transport, OGS resolves how water and solutes move through complex subsurface domains.
    visual:
      permalink: features/vis/DEMSelke3D.png
      alt: Groundwater model
  - eyebrow: 02 / Deep heat
    headline: Geothermal energy
    textline: Model heat extraction and storage in the deep subsurface, capturing the coupled thermal, hydraulic and mechanical response of geothermal reservoirs.
    visual:
      permalink: features/vis/layeredview2.png
      alt: Geothermal reservoir model
  - eyebrow: 03 / Long-term safety
    headline: Nuclear-waste disposal
    textline: Assess deep geological repositories over long timescales, resolving the full THMC evolution of the host rock and engineered barrier system.
    visual:
      permalink: features/vis/contours2-bw.png
      alt: Repository simulation contours

feature_intro: An adaptable, modular architecture enabling a wide variety of use cases and flexible research workflows.

features:
- headline: Pre-processing tools
  summary: Convert data sets, build and analyse meshes, and parametrize models with boundary conditions and source terms.
  textline: |
    A wide range of helper tools exist to get your model up and running with OpenGeoSys.

    Convert your existing data sets into appropriate OGS data formats and structures.

    Create meshes approximating geometrically the domain of interest. Analyze mesh quality, cleanup the mesh or adding layers to it.

    Parametrize the model with material parameters, boundary conditions and source terms.
  visual:
    permalink: "docs/tools/meshing-submeshes/extract-surface/TopBottomSideSurface.png"
    alt: Extracted surfaces

- headline: Process coupling
  summary: Solve coupled systems monolithically or with the staggered scheme for thermo-hydraulic, hydro-mechanical and phase-field problems.
  textline: |
    A coupled system of equations can be either solved in a fully coupled way of the monolithic method, or in the sequential manner of the staggered scheme. The monolithic scheme is applied for all coupled processes, while the staggered scheme are available for the coupled processes of thermo-hydraulic, hydro-mechanical, and phase field mechanical problems.
  visual:
    permalink: features/StaggeredCouplingScheme.png
    alt: Staggered coupling scheme

- headline: Data integration
  summary: Assess, integrate and visualize data sets with the OGS Data Explorer to catch artefacts and inconsistencies early.
  textline: |
    Integrate and visualize data sets for OpenGeoSys by using the OpenGeoSys Data Explorer. It provides functionality to visually assess the data and see possible artefacts, inconsistencies between data sets or missing information.
  visual:
    permalink: features/vis/chaohu_paper_mesh.png
    alt: OpenGeoSys mesh visualisation

- headline: Visualize results
  summary: Standard VTK output integrates directly with ParaView and VR-enabled visualization for intuitive exploration.
  textline: |
    By using VTK data formats visualizing simulation result data sets becomes an easy task. The de-facto standard software for scientific visualizations [ParaView](https://www.paraview.org) can be used to explore and analyze complex data in a visual way.

    [Virtual reality enabled visualization](https://www.ufz.de/vislab) brings your data onto the large screen for intuitive exploration and assessment.
  visual:
    permalink: features/vis/layeredview2.png
    alt: Layered model result visualisation

- headline: High-performance computing
  summary: Domain-decomposition parallelism built on PETSc and MPI scales across a wide variety of HPC architectures.
  textline: |
    High performance computing (HPC) has became a necessity in the modelling of environmental and geotechnical problems for better characterization of the complexity of geo-systems as well as predicting their evolution in time. Parallel computing is the most efficient method in the high performance computing. In OGS, the parallalization of the finite element (FE) computation is based on the domain decomposition method (DDC).

    Decomposed global matrices and vectors are handled by PETSc and the system of linear equations are solved by the performant PETSc solver. PETSc builds upon the Message Passing Interface (MPI) suitable for a wide variety of parallel computing architectures.

    Parallelization is implemented for single processes as well problems with coupled processes which are using the same order of element for each process.
  visual:
    permalink: "features/HPC-DDC.png"
    alt: Domain decomposition for parallel processing

- headline: Transparent development
  summary: A community-driven, fully open workflow with automated CI testing and mentored code review for every contribution.
  textline: |
    OpenGeoSys is an open-source project developed by a community of researchers. We try to be
    open-minded and and make team decisions. We try to help users and developers as best as we can.

    We invite you to take part in this journey, shape the future of OpenGeoSys together and happily welcome any contribution.
  visual:
    permalink: "features/OGS-Software-Engineering-Small.png"
    alt: Dev workflow

workflow_intro: The tools above are not isolated — they chain into a single reproducible pipeline. Every stage reads and writes open, standard formats, so a model flows from raw field data to published results without ever leaving the OGS ecosystem.

workflow_steps:
  - number: "01"
    headline: Input data
    textline: Field measurements, geological models and geometry imported into OGS data formats.
  - number: "02"
    headline: Mesh & parametrize
    textline: Build and refine meshes, then assign material parameters, boundary conditions and source terms.
  - number: "03"
    headline: Coupled simulation
    textline: Solve THMC processes monolithically or staggered, in parallel on HPC through PETSc and MPI.
  - number: "04"
    headline: Results & visualization
    textline: Export VTK for ParaView and VR — analysis, verification and publication-ready figures.

get_started:
  headline: Ready to dive in?
  items:
    - headline: Using OpenGeoSys
      textline: Start by downloading a prebuilt package for your platform and running the benchmarks.
      link:
        text: Download
        url: /releases
    - headline: Developing OpenGeoSys
      textline: Obtain the source, configure your toolchain and build the application from the Developer Guide.
      link:
        text: Developer Guide
        url: /docs/devguide
    - headline: Join the community
      textline: Get in touch through the Discourse forum, GitLab or email — every contribution is welcome.
      link:
        text: Discussion forum
        url: https://discourse.opengeosys.org
---
