/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on October 9, 2024, 6:03 PM
 */

#include "CheckMPLPhasesForSinglePhaseFlow.h"

#include <range/v3/algorithm/any_of.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/transform.hpp>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"

namespace MaterialPropertyLib
{
void checkMPLPhasesForSinglePhaseFlow(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map)
{
    // Check all of the elements have a medium defined.
    ranges::for_each(mesh.getElements() | MeshLib::views::ids,
                     [&](auto const& element_id)
                     { media_map.checkElementHasMedium(element_id); });

    // Collect phases of all elements...
    auto all_phases =
        media_map.media() |
        ranges::views::transform([&](auto const& medium)
                                 { return &fluidPhase(*medium); }) |
        ranges::to_vector;

    assert(!all_phases.empty());

    // ... and check if any of the phases are different by name.
    if (ranges::any_of(all_phases,
                       [p0 = all_phases.front()](auto const* const p)
                       { return p->name != p0->name; }))
    {
        OGS_FATAL(
            "You are mixing liquid and gas phases in your model domain. OGS "
            "does not yet know how to handle this.");
    }
}
}  // namespace MaterialPropertyLib
