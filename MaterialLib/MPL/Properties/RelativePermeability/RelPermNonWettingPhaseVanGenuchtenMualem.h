/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 1, 2020, 2:15 PM
 */

#pragma once

#include <limits>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;

class RelPermNonWettingPhaseVanGenuchtenMualem final : public Property
{
public:
    /**
     * @param name     Name of the property,
     * @param S_L_r    Residual saturation of the wetting phase,
     *                  \f$ S^L_r \f$
     * @param S_n_r    Residual saturation of the non-wetting phase,
     *                  \f$ S^n_{r} \f$
     * @param m        Exponent, \f$ m \in [0,1]\f$
     * @param krel_min Minimum relative permeability,
     *                  \f$ k_{rel}^n_{\mbox{min}}\f$
     */
    RelPermNonWettingPhaseVanGenuchtenMualem(std::string name,
                                             const double S_L_r,
                                             const double S_n_r,
                                             const double m,
                                             const double krel_min);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'RelPermNonWettingPhaseVanGenuchtenMualem' is "
                "implemented on the 'media' scale only.");
        }
    }

    /// \return \f$k_{rel}^n \f$.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;

    /// \return \f$ \frac{\mathrm{d} k_{rel}^n}{\mathrm{d}S^L} \f$.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    const double S_L_r_;     ///< Residual saturation of wetting phase.
    const double S_L_max_;   ///< Maximum saturation of wetting phase.
    const double m_;         ///< Exponent \f$ m \f$.
    const double krel_min_;  ///< Minimum relative permeability.
    const double
        S_L_for_krel_min_;  ///< Liquid saturation that gives \c krel_min_.

    /**
     *
     * @param S_L Liquid saturation.
     * @return \return \\return \f$k_{rel}^n \f$.
     */
    double computeValue(const double S_L) const;
    /**
     *  Computes
     * \f[\frac{\mathrm{d} k_{rel}^n}{\mathrm{d}S^L}=
     * -(\dfrac{[1-S_e^{1/m}]^{2m}}{2\sqrt{1-S_e}}+2\sqrt{1-S_e}
     * {(1-S_e^{1/m})}^{2*m-1}
     *  S_e^{1/m-1})/(S^L_\mbox{max}-S^L_r) \f]
     *  As \f$S^L \to S^L_\mbox{max}\f$, or \f$S_e \to 1\f$,
     * \f$\dfrac{[1-S_e^{1/m}]^{2m}}{2\sqrt{1-S_e}}\f$ has
     * limit zero.
     * @param S_L Liquid saturation.
     * @return \return \f$ \frac{\mathrm{d} k_{rel}^n}{\mathrm{d}S^L} \f$.
     */
    double computeDerivative(const double S_L) const;
    /**
     * Computes the saturation that gives the minimum relative permeability by
     * using the Newton-Raphson method.
     * @return \f$ S^L\f$ that gives the minimum relative permeability.
     */
    double computeSaturationForMinimumRelativePermeability() const;
};
}  // namespace MaterialPropertyLib
