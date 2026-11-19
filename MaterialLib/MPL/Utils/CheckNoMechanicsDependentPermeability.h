// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string_view>

namespace MaterialPropertyLib
{
class MaterialSpatialDistributionMap;

/// Refuses a permeability that depends on the stress or on the strain in a
/// process that does not solve the deformation equation, by calling OGS_FATAL,
/// which throws std::runtime_error. Media without a permeability are skipped.
///
/// The refused properties are
///   - StrainDependentPermeability,
///   - PermeabilityMohrCoulombFailureIndexModel,
///   - EmbeddedFracturePermeability,
///   - OrthotropicEmbeddedFracturePermeability.
///
/// \param media_map    the media of the process.
/// \param process_name the process type, used in the error message only.
void checkNoMechanicsDependentPermeability(
    MaterialSpatialDistributionMap const& media_map,
    std::string_view process_name);
}  // namespace MaterialPropertyLib
