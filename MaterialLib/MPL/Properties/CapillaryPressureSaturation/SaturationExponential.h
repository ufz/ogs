/**
 * \file
 * \author Norbert Grunwald
 * \date   27.06.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
/**
 * \brief A simplistic soil characteristics function.
 * \details This property must be a medium property. It
 * computes the saturation of the wetting phase as function
 * of capillary pressure.
 *
 * The capillary pressure-saturation relationship given by
 * \f$s_{eff}=1-\left(\frac{p_{c}}{p_{c}^{ref}}\right)^{\lambda}\f$,
 * where
 * \f$\lambda\f$ is an exponent,
 * \f$p_{c}\f$ is capillary pressure,
 * \f$p_{c}^{ref}\f$ is reference capillary pressure at which \f$s_{eff}=0\f$.
 *
 * This property can mainly be used for testing. If the exponent is set to 1,
 * the characteristic curve shows a linear dependence.
 *
 *
 * Reference capillary pressure at which \f$s_{eff}=0\f$.
 */
class SaturationExponential final : public Property
{
private:
    const double
        residual_liquid_saturation_;  ///< Residual saturation of the gas phase.
    const double
        residual_gas_saturation_;  ///< Residual saturation of the liquid phase.
    const double p_cap_ref_;       ///< Reference capillary pressure at which
                                   ///  effective saturation is zero.
    const double exponent_;        ///< Exponent to govern the shape of the curve.

public:
    SaturationExponential(std::string name,
                          const double residual_liquid_saturation,
                          const double residual_gas_saturation,
                          const double p_cap_ref, const double exponent);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'SaturationExponential' is implemented on the "
                "'media' scale only.");
        }
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& /*pos*/,
                           double const /*t*/,
                           double const /*dt*/) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override;
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const /*primary_variable1*/,
                             Variable const /*primary_variable2*/,
                             ParameterLib::SpatialPosition const& /*pos*/,
                             double const /*t*/,
                             double const /*dt*/) const override;
};
}  // namespace MaterialPropertyLib
