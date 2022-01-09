/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PhaseTransitionModel.h"

namespace ProcessLib
{
namespace TH2M
{
int numberOfComponents(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media,
    std::string phase_name)
{
    // Always the first (begin) medium that holds fluid phases.
    auto const medium = media.begin()->second;
    return medium->phase(phase_name).numberOfComponents();
}

int findComponentIndex(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media,
    std::string phase_name, MaterialPropertyLib::PropertyType property_type)
{
    // It is always the first (begin) medium that holds fluid phases.
    auto const medium = media.begin()->second;
    auto const& phase = medium->phase(phase_name);

    // find the component for which the property 'property_type' is defined
    for (std::size_t c = 0; c < phase.numberOfComponents(); c++)
    {
        if (phase.component(c).hasProperty(property_type))
        {
            return c;
        }
    }

    // A lot of checks can (and should) be done to make sure that the right
    // components with the right properties are used. For example, the names
    // of the components can be compared to check that the name of the
    // evaporable component does not also correspond to the name of the
    // solvate.

    OGS_FATAL(
        "PhaseTransitionModel: findComponentIndex() could not find the "
        "specified property type '{:s}' in phase '{:s}'.",
        MaterialPropertyLib::property_enum_to_string[property_type],
        phase_name);
}
}  // namespace TH2M
}  // namespace ProcessLib