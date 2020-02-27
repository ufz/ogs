/**
 * \file
 * \author Norbert Grunwald
 * \date   27.06.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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

/// Orthotropic, saturation dependent swelling model.
/// This property must be a solid phase property, it computes the stress
/// increment depending on the saturation.  A local coordinate system can be
/// given for orthotropy.
class SaturationDependentSwelling final : public Property
{
private:
    Phase* _phase = nullptr;
    /// Maximum swelling pressures, one for each spatial dimension.
    std::array<double, 3> const _p;
    /// Exponents, one for each spatial dimension.
    std::array<double, 3> const _lambda;
    double const _S_min;  //< Lower saturation limit.
    double const _S_max;  //< Upper saturation limit.
    ParameterLib::CoordinateSystem const* const _local_coordinate_system;

public:
    SaturationDependentSwelling(
        std::array<double, 3> swelling_pressures,
        std::array<double, 3>
            exponents,
        double const lower_saturation_limit,
        double const upper_saturation_limit,
        ParameterLib::CoordinateSystem const* const local_coordinate_system);

    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override;

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};

}  // namespace MaterialPropertyLib
