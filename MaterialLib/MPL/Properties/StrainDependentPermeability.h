/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 10, 2020, 8:49 AM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"

namespace ParameterLib
{
struct CoordinateSystem;
template <typename T>
struct Parameter;
}  // namespace ParameterLib

namespace MaterialPropertyLib
{
template <int DisplacementDim>
class StrainDependentPermeability final : public Property
{
public:
    StrainDependentPermeability(
        std::string name, ParameterLib::Parameter<double> const& k0,
        double const b1, double const b2, double const b3,
        double const minimum_permeability, double const maximum_permeability,
        ParameterLib::CoordinateSystem const* const local_coordinate_system);

    void checkScale() const override;

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    /// Initial intrinsic permeability.
    ParameterLib::Parameter<double> const& k0_;
    double const b1_;
    double const b2_;
    double const b3_;
    double const minimum_permeability_;
    double const maximum_permeability_;
    ParameterLib::CoordinateSystem const* const local_coordinate_system_;
};

extern template class StrainDependentPermeability<2>;
extern template class StrainDependentPermeability<3>;

}  // namespace MaterialPropertyLib
