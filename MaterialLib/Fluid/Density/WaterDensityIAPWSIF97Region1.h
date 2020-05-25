/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on December 8, 2016, 4:19 PM
 */

#pragma once

#include <string>

#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/GibbsFreeEnergy/DimensionLessGibbsFreeEnergyRegion1.h"

namespace MaterialLib
{
namespace Fluid
{
/** Water density model
 *  base on the IAPWS Industrial Formulation 1997
 *  <a href="http://www.iapws.org/relguide/IF97-Rev.pdf">IF97-Rev</a>
 */
class WaterDensityIAPWSIF97Region1 final : public FluidProperty
{
public:
    WaterDensityIAPWSIF97Region1()
        : gibbs_free_energy_(DimensionLessGibbsFreeEnergyRegion1()){};

    /// Get density model name.
    std::string getName() const override
    {
        return "Water density model based on IAPWSIF97";
    }

    /// Get density value.
    /// \param var_vals Variable values in an array of temperature and pressure.
    double getValue(const ArrayType& var_vals) const override
    {
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        const double p = var_vals[static_cast<int>(PropertyVariableType::p)];
        const double tau = ref_T_ / T;
        const double pi = p / ref_p_;

        return ref_p_ / (sR_ * T * gibbs_free_energy_.get_dgamma_dpi(tau, pi));
    }

    /**
     *  Get the partial differential of the density with respect to temperature
     *  or pressure.
     *  \param var_vals Variable values in an array of temperature and pressure.
     *  \param var_type Variable type.
     */
    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var_type) const override
    {
        const double T = var_vals[static_cast<int>(PropertyVariableType::T)];
        const double p = var_vals[static_cast<int>(PropertyVariableType::p)];

        const double tau = ref_T_ / T;
        const double pi = p / ref_p_;

        const double dgamma_dpi = gibbs_free_energy_.get_dgamma_dpi(tau, pi);
        switch (var_type)
        {
            case PropertyVariableType::T:
                return -(ref_p_ -
                         tau * ref_p_ *
                             gibbs_free_energy_.get_dgamma_dtau_dpi(tau, pi) /
                             dgamma_dpi) /
                       (sR_ * T * T * dgamma_dpi);
            case PropertyVariableType::p:
                return -gibbs_free_energy_.get_dgamma_dpi_dpi(tau, pi) /
                       (sR_ * T * dgamma_dpi * dgamma_dpi);
            default:
                return 0.0;
        }
    }

private:
    const DimensionLessGibbsFreeEnergyRegion1 gibbs_free_energy_;

    const double ref_T_ = 1386;     ///< reference temperature in K.
    const double ref_p_ = 1.653e7;  ///< reference pressure in Pa.
    /// Specific water vapour gas constant in J/(kgK).
    const double sR_ = 461.526;
};

}  // namespace Fluid
}  // namespace MaterialLib
