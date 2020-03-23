/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
class Phase;
class Component;

/**
 *  \brief This class handles the computation of \f$ p_c(S) \f$ with the
 * capillary pressure regularization.
 *
 *   For the regularized van Genuchten model, please refer to
 *  \cite marchand2013fully.
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

    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override
    {
        if (!std::holds_alternative<Medium*>(scale_pointer))
        {
            OGS_FATAL(
                "The property 'CapillaryPressureRegularizedVanGenuchten' is "
                "implemented on the "
                "'media' scale only.");
        }
        _medium = std::get<Medium*>(scale_pointer);
    }

    /// It returns \f$ p_c(S) \f$.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    /// It returns \f$ \frac{\partial p_c(S)}{\partial  S} \f$
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;

    /// It returns \f$ \frac{\partial^2 p_c(S)}{\partial  S^2} \f$
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1, Variable const variable2,
                             ParameterLib::SpatialPosition const& pos,
                             double const t, double const dt) const override;

private:
    Medium* _medium = nullptr;
    const double _saturation_nonwet_r;  ///< Residual saturation of nonwetting
    /// Residual saturation of liquid phase.
    double const _residual_saturation;
    double const _maximuml_saturation;  ///< Maximum saturation of liquid phase.
    double const _m;                    ///< Exponent.
    /// Capillary pressure scaling factor. Sometimes, it is called apparent gas
    /// entry pressure.
    double const _p_b;
    /// parameter in regularized van Genuchten model
    static constexpr double _xi = 1e-5;

    /// Gets regularized capillary pressure via the saturation of the
    /// non-wetting phase.
    double getPcBarvGSg(double const Sg) const;
    /// Gets regularized saturation via the saturation of the non-wetting
    /// phase.
    double getSBar(double const Sg) const;
    /// Gets capillary pressure via the non-wetting phase saturation with the
    /// original van Genuchten model.
    double getPcvGSg(double const Sg) const;
    /// Gets \f$\dfrac{\partial}{\partial {\bar S}_g }\f$.
    double getdPcdSvGBar(double const Sg) const;
    /// Gets \f$\dfrac{\partial}{\partial {S}_g }\f$.
    double getdPcdSvG(double const Sg) const;
};

}  // namespace MaterialPropertyLib
