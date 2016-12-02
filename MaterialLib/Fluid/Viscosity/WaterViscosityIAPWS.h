/**
 *  \brief Viscosity model according to
 *         Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary
 *         Water Substance
 *
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   WaterViscosityIAPWS.h
 *
 * Created on December 1, 2016, 1:41 PM
 */

#ifndef OGS_WATER_VISCOSITY_IAPWS_H
#define OGS_WATER_VISCOSITY_IAPWS_H

#include <cmath>
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

    /** Get the partial differential of the density with respect to temperatur
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

    const double _hi[4] = {1.67752, 2.20462, 0.6366564, -0.241605};
    const double _hij[6][7] = {
        {0.520094, 0.222531, -0.281378, 0.161913, -0.0325372, 0, 0},
        {0.0850895, 0.999115, -0.906851, 0.257399, 0, 0, 0},
        {-1.08374, 1.88797, -0.772479, 0, 0, 0, 0},
        {-0.289555, 1.26613, -0.489837, 0, 0.0698452, 0, -0.00435673},
        {0, 0, -0.25704, 0, 0, 0.00872102, 0},
        {0, 0.120573, 0, 0, 0, 0, -0.000593264}};

    double computeBarMu0Factor(const double barT) const;
    double computeBarMu1Factor(const double barT, const double bar_rho) const;

    double computedBarMu_dbarT(const double barT, double bar_rho) const;
    double computedBarMu_dbarRho(const double barT, double bar_rho) const;
};

}  // end namespace
}  // end namespace

#endif /* OGS_WATER_VISCOSITY_IAPWS_H */
