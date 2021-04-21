/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "MaterialLib/MPL/Property.h"

namespace ParameterLib
{
struct CoordinateSystem;
template <typename T>
struct Parameter;
}  // namespace ParameterLib

namespace MaterialPropertyLib
{
class Medium;
/// Porosity mixing based model for effective heat conduction
/// \details This property is a medium property.
/// The corresponding values are taken from the liquid/solid phase.
template <int GlobalDim>
class EffectiveThermalConductivityPorosityMixing final : public Property
{
public:
    EffectiveThermalConductivityPorosityMixing(std::string name,
            ParameterLib::CoordinateSystem const* const local_coordinate_system);

    void checkScale() const override;

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override;

private:
    ParameterLib::CoordinateSystem const* const local_coordinate_system_;
};

extern template class EffectiveThermalConductivityPorosityMixing<2>;
extern template class EffectiveThermalConductivityPorosityMixing<3>;

}  // namespace MaterialPropertyLib
