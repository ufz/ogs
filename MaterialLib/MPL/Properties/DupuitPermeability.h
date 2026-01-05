// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"
#include "ParameterLib/Parameter.h"

namespace MaterialPropertyLib
{
class DupuitPermeability final : public Property
{
public:
    /// This constructor accepts two parameters.
    DupuitPermeability(std::string name,
                       ParameterLib::Parameter<double> const& parameter);

    /// This method computes the value of a property depending linearly on
    /// the value of the given primary variable.
    PropertyDataType value(
        MaterialPropertyLib::VariableArray const& variable_array,
        ParameterLib::SpatialPosition const& pos, double const t,
        double const dt) const override;

private:
    ParameterLib::Parameter<double> const& parameter_;
};
}  // namespace MaterialPropertyLib
