/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   WaterDensityIAPWSIF97Region1.h
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
        : _gibbs_free_energy(DimensionLessGibbsFreeEnergyRegion1()){};

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
        const double tau = _ref_T / T;
        const double pi = p / _ref_p;

        return _ref_p / (_sR * T * _gibbs_free_energy.get_dgamma_dpi(tau, pi));
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

        const double tau = _ref_T / T;
        const double pi = p / _ref_p;

        const double dgamma_dpi = _gibbs_free_energy.get_dgamma_dpi(tau, pi);
        switch (var_type)
        {
            case PropertyVariableType::T:
                return -(_ref_p -
                         tau * _ref_p *
                             _gibbs_free_energy.get_dgamma_dtau_dpi(tau, pi) /
                             dgamma_dpi) /
                       (_sR * T * T * dgamma_dpi);
            case PropertyVariableType::p:
                return -_gibbs_free_energy.get_dgamma_dpi_dpi(tau, pi) /
                       (_sR * T * dgamma_dpi * dgamma_dpi);
            default:
                return 0.0;
        }
    }

private:
    const DimensionLessGibbsFreeEnergyRegion1 _gibbs_free_energy;

    const double _ref_T = 1386;     ///< reference temperature in K.
    const double _ref_p = 1.653e7;  ///< reference pressure in Pa.
    /// Specific water vapour gas constant in J/(kgK).
    const double _sR = 461.526;
};

}  // end of namespace
}  // end of namespace
