/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 4, 2021, 3:05 PM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Phase;

/**
 *  \brief A model for water vapour density in the unsaturated porous media.
 *
 *   The water vapour density is given by
 *   \f[
 *      \rho_v = h\, \rho_{vS},
 *   \f]
 *  where \f$h\f$ is the relative humidity according to
 *   \f[
 *     h=\exp({\frac{p}{\rho_w R T}}),
 *   \f]
 *
 *  with \f$R\f$ the specific gas constant for water vapour, and
 * \f$\rho_{vS}\f$, is the saturated vapour density given by
 *   \f[
 *     \rho_{vS}=10^{-3}\, \exp({19.819-4975.9/T}).
 *   \f]
 * The specific gas constant for water vapour
 * \f$R=461.6\,\text{J}\,\text{kg}^{-1}\text{K}^{-1}\f$,
 * the value is calculated by \f$R_g/M_w\f$ with
 * \f$R_g= 8.3144621\,\text{J}\,\text{mol}^{-1}\text{K}^{-1}\f$,
 *  the ideal gas constant, and
 * \f$M_w=0.018016\,\text{kg}\,\text{mol}^{-1}\f$,
 *  the molar mass of water.
 *
 *  The formula of this density model is presented on page 20 of
 *  \cite kimball1976comparison.
 *
 *  In the application, the vapour density related terms in the mass balance
 *   equation are multiplied with \f$1-S\f$ with \f$ S \f$, the water
 * saturation. Therefore the application of the water vapour density model is
 * naturally restricted in the unsaturated zones.
 */
class WaterVapourDensity final : public Property
{
public:
    explicit WaterVapourDensity(std::string name) { name_ = std::move(name); }
    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'WaterVapourDensity' is "
                "implemented on the 'Phase' scale only.");
        }
    }

    /// \return The water vapour density.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    /// \return The derivative of  water vapour density with respect to
    /// temperature or phase (water) pressure.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};
}  // namespace MaterialPropertyLib
