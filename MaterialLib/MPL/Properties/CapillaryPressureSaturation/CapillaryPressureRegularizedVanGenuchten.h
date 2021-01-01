/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 20, 2020, 9:30 AM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;

/**
 *  \brief This class handles the computation of the capillary pressure,
 *  \f$ p_c(S_l) \f$, with the
 * capillary pressure regularization.
 *
 *   For the regularized van Genuchten model, please refer to
 *  \cite huang2015extending.
 *
 *   For the van Genuchten capillary pressure mode,
 *   \sa MaterialPropertyLib::CapillaryPressureVanGenuchten
 */
class CapillaryPressureRegularizedVanGenuchten final : public Property
{
public:
    CapillaryPressureRegularizedVanGenuchten(
        double const residual_liquid_saturation,
        double const maximum_liquid_saturation,
        double const exponent,
        double const p_b);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'CapillaryPressureRegularizedVanGenuchten' is "
                "implemented on the 'media' scale only.");
        }
    }

    /// \return \f$ p_c(S_l) \f$.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    /// \return \f$ \frac{\partial p_c(S_l)}{\partial  S_l} \f$
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;

private:
    double const Sg_r_;    ///< Residual saturation of gas phase
    double const Sg_max_;  ///< Maximum saturation of gas phase
    double const m_;       ///< Exponent.
    /// Capillary pressure scaling factor. Sometimes, it is called apparent gas
    /// entry pressure.
    double const p_b_;
    /// parameter in regularized van Genuchten model
    static constexpr double xi_ = 1e-5;

    double const PcBarvGSg_Sg_max_;
    double const dPcdSvGBarSg_max_;

    /// Gets regularized capillary pressure via the saturation of the
    /// non-wetting phase.
    double getPcBarvGSg(double const Sg) const;
    /// Gets regularized saturation via the saturation of the non-wetting
    /// phase.
    double getSBar(double const Sg) const;
    /// Gets capillary pressure via the non-wetting phase saturation with the
    /// original van Genuchten model.
    double getPcvGSg(double const Sg) const;
    /// Gets \f$\frac{\partial p_c}{\partial {\bar S}_g }\f$.
    double getdPcdSvGBar(double const Sg) const;
    /// Gets \f$\frac{\partial p_c}{\partial {S}_g }\f$.
    double getdPcdSvG(double const Sg) const;
};

}  // namespace MaterialPropertyLib
