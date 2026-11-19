// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CheckNoMechanicsDependentPermeability.h"

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/EmbeddedFracturePermeability.h"
#include "MaterialLib/MPL/Properties/OrthotropicEmbeddedFracturePermeability.h"
#include "MaterialLib/MPL/Properties/PermeabilityMohrCoulombFailureIndexModel.h"
#include "MaterialLib/MPL/Properties/StrainDependentPermeability.h"

namespace MaterialPropertyLib
{
namespace
{
/// The properties exist in 2D and 3D only.
template <template <int> class PropertyClass>
bool isOfType(Property const& property)
{
    return dynamic_cast<PropertyClass<2> const*>(&property) ||
           dynamic_cast<PropertyClass<3> const*>(&property);
}

/// Returns the type name of a permeability that reads the stress or the
/// strain, or nullptr if the permeability reads neither.
char const* mechanicsDependentTypeName(Property const& permeability)
{
    if (isOfType<StrainDependentPermeability>(permeability))
    {
        return "StrainDependentPermeability";
    }
    if (isOfType<PermeabilityMohrCoulombFailureIndexModel>(permeability))
    {
        return "PermeabilityMohrCoulombFailureIndexModel";
    }
    if (isOfType<EmbeddedFracturePermeability>(permeability))
    {
        return "EmbeddedFracturePermeability";
    }
    if (isOfType<OrthotropicEmbeddedFracturePermeability>(permeability))
    {
        return "OrthotropicEmbeddedFracturePermeability";
    }
    return nullptr;
}
}  // namespace

void checkNoMechanicsDependentPermeability(
    MaterialSpatialDistributionMap const& media_map,
    std::string_view const process_name)
{
    for (auto const& medium : media_map.media())
    {
        // A medium that no element uses need not define a permeability.
        if (!medium->hasProperty(PropertyType::permeability))
        {
            continue;
        }
        if (auto const* const type_name = mechanicsDependentTypeName(
                medium->property(PropertyType::permeability)))
        {
            OGS_FATAL(
                "The permeability of {:s} is of type '{:s}', which depends on "
                "the stress or the strain. The {:s} process does not solve the "
                "deformation equation and computes neither. Use a "
                "permeability that does not depend on the mechanics.",
                medium->description(), type_name, process_name);
        }
    }
}
}  // namespace MaterialPropertyLib
