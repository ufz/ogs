// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    EffectiveThermalConductivityPorosityMixing(
        std::string name,
        ParameterLib::CoordinateSystem const* const local_coordinate_system);

    void checkScale() const override;

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;

private:
    ParameterLib::CoordinateSystem const* const local_coordinate_system_;
};

extern template class EffectiveThermalConductivityPorosityMixing<1>;
extern template class EffectiveThermalConductivityPorosityMixing<2>;
extern template class EffectiveThermalConductivityPorosityMixing<3>;

}  // namespace MaterialPropertyLib
