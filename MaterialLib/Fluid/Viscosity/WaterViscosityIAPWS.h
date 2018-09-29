/**
 *  \brief Viscosity model according to
 *         Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary
 *         Water Substance
 *
 *  \copyright
 *   Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   WaterViscosityIAPWS.h
 *
 * Created on December 1, 2016, 1:41 PM
 */

#pragma once

#include <string>

#include "MaterialLib/Fluid/FluidProperty.h"

namespace MaterialLib
{
namespace Fluid
{
/**
 * \brief A class for viscosity model that is defined by
 *        The International Association for the Properties of Water and Steam
 *        <a href="http://www.iapws.org/relguide/visc.pdf">IAPWS</a>
 *
 *        With the definition, the viscosity is a function of temperature and
 *        water density
 *
 *  \attention The critical enhancement, \f$\bar{\mu}_2\f$, which is significant
 *             only within the boundaries specified by
 *                 \f[ T (\mbox{in K}) \in (645.91, 650.77) \f]
 *             and
 *                 \f[ \rho (\mbox{in kg m}^{-3}) \in (245.8, 405.3)\f],
 *             is not considered.
 */
class WaterViscosityIAPWS final : public FluidProperty
{
public:
    WaterViscosityIAPWS() = default;

    /// Get model name.
    std::string getName() const override
    {
        return "IAPWS temperature-density dependent viscosity model";
    }

    /**
     *  Get density value.
     *  \param var_vals Variable values of temperature and water density in
     *                  an array. The order of its elements
     *                  is given in enum class PropertyVariableType.
     *
     *  \attention The variable values must follow the SI standard,
     *             e.g temperature unit is K.
     */
    double getValue(const ArrayType& var_vals) const override;

    /** Get the partial differential of the density with respect to temperature
     *  or liquid pressure.
     *  \param var_vals Variable values of temperature and water density in
     *                  an array. The order of its elements
     *                  is given in enum class PropertyVariableType.
     *
     *  \attention The variable values must follow the SI standard,
     *             e.g temperature unit is K.
     *  \param var_type Variable type.
     */
    double getdValue(const ArrayType& var_vals,
                     const PropertyVariableType var_type) const override;

private:
    const double _ref_T = 647.096;  ///< reference temperature in K
    const double _ref_rho = 322.0;  ///< reference density in `kg/m^3`
    const double _ref_mu = 1.0e-6;  ///< reference viscosity in Pa.s

    // Coefficients Hi and Hij are given in two static arrays in the cpp file.
};

}  // end namespace
}  // end namespace
