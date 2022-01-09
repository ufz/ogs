/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 19, 2021, 11:49 AM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
class Phase;

/**
 * \brief A latent heat model of vaporisation of water considering the critical
 *  temperature.
 *
 *  The model uses an equation for a general expression of the latent heat of
 *  vaporisation of water in the vicinity of and far away from the critical
 *  temperature, which was presented by Torquato and Stell in
 *  \cite torquato1982equation.
 *
 *  Denoting the critical temperature as \f$T_c\f$, and introducing a
 *  dimensionless variable \f$\tau=(T_c-T)/T_c\f$ associated with temperature
 *  \f$T\f$, the equation is given by
 *
 *  \f[
 *    L(\tau) = a_1 \tau^{\beta}+a_2 \tau^{\beta+\Delta}
 *              +a_4 \tau^{1-\alpha+\beta}
 *             +\sum_{n=1}^{M}(b_n \tau^n),\,\text{[kJ/kg]},
 *  \f]
 *   where the parameters of \f$b_n\f$ are obtained by the least square method
 *   by fitting the equation with the experiment data.
 *
 *   In this model, the parameter set of \f$M=5\f$ is taken for a high accuracy.
 *   All parameters are given below:
 *  <ul>
 *    <li> \f$\alpha=1/8,\,\beta=1/3,\, \Delta=0.79-\beta\f$,
 *    <li> \f$a_1=1989.41582,\, a_2=11178.45586, a_4=26923.68994\f$,
 *    <li> \f$b_n:=\{-28989.28947, -19797.03646, 28403.32283,
 *                                  -30382.306422, 15210.380\}\f$.
 *  </ul>
 *
 *  The critical temperature is 373.92 \f$^{\circ}\f$C.
 *
 *  A comparison of this model with the model of
 *  MaterialPropertyLib::LinearWaterVapourLatentHeat
 *  is given in the following figure.
 *
 *  \image{inline} html latent_heat_water_vapour_with_critical_T.png ""
 *
 */
class WaterVapourLatentHeatWithCriticalTemperature final : public Property
{
public:
    explicit WaterVapourLatentHeatWithCriticalTemperature(std::string name)
    {
        name_ = std::move(name);
    }

    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'WaterVapourLatentHeatWithCriticalTemperature' "
                "is "
                "implemented on the 'phase' scale only.");
        }
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};

}  // namespace MaterialPropertyLib
