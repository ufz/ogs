/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 20, 2022
 */
#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/SigmoidFunction.h"
namespace MaterialPropertyLib
{
/**
 * Temperature dependent model for some volume or mass fraction
 *
 * \details This property must be a medium property.
 * It can be used for media with a phase change of the fluid
 * in the pore space. Then, there is a temperature
 * dependent (volume or mass) fraction for one of the two phases.
 * The other one is immediately given by the remaining fraction.
 * The fraction is computed based on a phase transition spread
 * over a temperature interval following a sigmoid line:
 *
 * \f[
 *      \frac{\phi}{1 + \exp(k(T - T_\mathrm{c}))}
 * \f]
 *
 * Parameter \f$k\f$ controlling the steepness of the phase change
 * and \f$T_\mathrm{c}\f$ as the characteristic temperature where
 * the phase change is centered around (melting temperature)
 * are forwarded to a sigmoid function.
 *
 **/
class TemperatureDependentFraction final : public Property
{
public:
    TemperatureDependentFraction(std::string name,
                                 double const k,
                                 double const T_c);

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

    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1, Variable const variable2,
                             ParameterLib::SpatialPosition const& pos,
                             double const t, double const dt) const override;

private:
    SigmoidFunction phase_change_shape_;  //< shape of the phase change
};

}  // namespace MaterialPropertyLib
