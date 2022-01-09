/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
/// Bishop's effective stress model using saturation cutoff. The effective
/// stress is 1 as long as the saturation does not fall below the saturation
/// cutoff value.
class BishopsSaturationCutoff final : public Property
{
public:
    BishopsSaturationCutoff(std::string name, double const cutoff_value);

    void checkScale() const override;

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& /*pos*/,
                           double const /*t*/,
                           double const /*dt*/) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override;

private:
    double const S_L_max_;  //< Maximum saturation cutoff value.
};
}  // namespace MaterialPropertyLib
