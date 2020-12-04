/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
/// saturation dependend model for effective heat conduction
/// \details This property must be a medium property, it
/// computes the effetive heat conductivity based on a wet
/// and a dry value
class HeatConductionSaturationDependent final : public Property
{
public:
    HeatConductionSaturationDependent(std::string name,
                                      double const K_dry,
                                      double const K_wet);

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
    double const K_dry_;
    double const K_wet_;
};
}  // namespace MaterialPropertyLib
