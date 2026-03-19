// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "BoreholeGeometry.h"

#include "BHESectionUtils.h"
#include "BaseLib/ConfigTree.h"
#include "MeshLib/Node.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/SpatialPosition.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
static BoreholeGeometry buildBoreholeGeometry(
    double borehole_length,
    ParameterLib::Parameter<double> const& diameter_parameter,
    std::vector<MeshLib::Node*> const& bhe_nodes)
{
    auto const distances = cumulativeDistances(bhe_nodes);

    // Sample diameters at each node
    std::vector<double> diameters;
    diameters.reserve(bhe_nodes.size());
    for (auto const* const node : bhe_nodes)
    {
        ParameterLib::SpatialPosition pos;
        pos.setNodeID(node->getID());
        pos.setCoordinates(*node);
        auto const diameter = diameter_parameter(0.0 /*t*/, pos)[0];
        if (diameter <= 0)
        {
            OGS_FATAL(
                "Non-positive borehole diameter {:g} at node {:d} "
                "(x={:g}, y={:g}, z={:g}).",
                diameter, node->getID(), (*node)[0], (*node)[1], (*node)[2]);
        }
        diameters.push_back(diameter);
    }

    auto [section_boundaries, section_diameters] =
        groupSections(distances, diameters);

    return {borehole_length,
            {std::move(section_boundaries), std::move(section_diameters)},
            &diameter_parameter};
}

BoreholeGeometry BoreholeGeometry::rebuildForNodes(
    std::vector<MeshLib::Node*> const& new_nodes) const
{
    if (!diameter_param)
    {
        OGS_FATAL(
            "Cannot rebuild borehole geometry: no diameter parameter stored. "
            "This geometry was not created from a parameter-based config.");
    }
    return buildBoreholeGeometry(length, *diameter_param, new_nodes);
}

BoreholeGeometry createBoreholeGeometry(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::vector<MeshLib::Node*> const& bhe_nodes)
{
    const auto borehole_length =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__borehole__length}
        config.getConfigParameter<double>("length");
    if (borehole_length <= 0)
    {
        OGS_FATAL("Borehole length must be positive, got {:g}.",
                  borehole_length);
    }

    auto const diameter_parameter_or_value =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__borehole__diameter}
        config.getConfigParameter<std::string>("diameter");
    auto const& diameter_parameter =
        ParameterLib::getNamedOrCreateInlineParameter(
            diameter_parameter_or_value, parameters, "borehole_diameter",
            "borehole_diameter");

    return buildBoreholeGeometry(borehole_length, diameter_parameter,
                                 bhe_nodes);
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
