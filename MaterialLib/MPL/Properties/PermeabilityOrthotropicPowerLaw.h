/**
 * \file
 * \author Norbert Grunwald
 * \date   27.06.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "MaterialLib/MPL/Property.h"


namespace ParameterLib
{
struct CoordinateSystem;
}

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

/// Orthotropic permeability model based on power law dependency on porosity.
/// \details This property must be a solid phase property, it
/// computes the permeability depending on the porosity. A local coordinate
/// system can be given for orthotropy.
template <int DisplacementDim>
class PermeabilityOrthotropicPowerLaw final : public Property
{
private:
    /// Intrinsic permeabilities, one for each spatial dimension.
    std::array<double, DisplacementDim> const k_;
    /// Exponents, one for each spatial dimension.
    std::array<double, DisplacementDim> const lambda_;
    ParameterLib::CoordinateSystem const* const local_coordinate_system_;

public:
    PermeabilityOrthotropicPowerLaw(
        std::string name,
        std::array<double, DisplacementDim>
            intrinsic_permeabilities,
        std::array<double, DisplacementDim>
            exponents,
        ParameterLib::CoordinateSystem const* const local_coordinate_system);

    void checkScale() const override;

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};

}  // namespace MaterialPropertyLib
