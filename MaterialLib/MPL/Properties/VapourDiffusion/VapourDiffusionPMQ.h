/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 7, 2021, 9:14 AM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Phase;

/**
 * \brief The Penman-Millington-Quirk (PMQ) Vapour diffusion model.
 *
 *  The vapour diffusion can be described by
 *   \cite moldrup1997modeling, \cite moldrup2000predicting,
 *  \f[
 * 	D_v=D_0 \left(\frac{T}{273.15}\right)^{n}
 * 	D_{vr},
 *  \f]
 *  where \f$D_{0}\f$ is the base diffusion coefficient with default value
 *  \f$2.16\cdot 10^{-5}\f$ \f${\text m}^2
 * \text{Pa}/(\text{s}\text{K}^{n})\f$,
 *  \f$n\f$ is the exponent with default value 1.8,
 *  \f$D_{vr}\f$ is the the relative
 * diffusion coefficient, and \f$T\f$ is the temperature.
 *
 *  The Penman–Millington–Quirk (PMQ) model \cite moldrup1997modeling is given
 *  as
 *  \f[
 *     D_{vr}=0.66 \phi \left(\frac{\kappa}{\phi}\right)^{\frac{12-m}{3}},
 *  \f]
 * where \f$\phi\f$ is the total porosity, \f$\kappa\f$ is the air filled
 *  porosity, and \f$m\f$ is a fitting parameter. The air filled porosity is
 *  defined as \f$ \kappa = \phi-\theta = \phi -S_L \phi \f$ with \f$\theta\f$
 *  the liquid content, and \f$S_L\f$ the liquid saturation.
 *
 * According to the study presented in \cite moldrup1997modeling, \f$m=6\f$
 *  is the best fitting parameter for the sieved, repacked soils that the
 *  authors tested. Therefore,  \f$m=6\f$ is used in the implementation, which
 * gives
 *   \f[
 *     D_{vr}=0.66 \phi (1 - S_L )^2.
 *   \f]
 *
 * Note: In order to maintain consistency with the implementation of the
 *  computations of other vapor-related parameters, \f$ \phi (1 - S_L )\f$ is
 *  removed from the implementation for this class and is multiplied back in
 *  the local assembler.
 */
class VapourDiffusionPMQ final : public Property
{
public:
    explicit VapourDiffusionPMQ(std::string name,
                                double const base_diffusion_coefficient,
                                double const exponent)
        : base_diffusion_coefficient_(base_diffusion_coefficient),
          exponent_(exponent)
    {
        name_ = std::move(name);
    }

    void checkScale() const override
    {
        if (!(std::holds_alternative<Phase*>(scale_) ||
              std::holds_alternative<Component*>(scale_)))
        {
            OGS_FATAL(
                "The property 'VapourDiffusionPMQ' is "
                "implemented on the 'phase' and 'component' scale only.");
        }
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    double const base_diffusion_coefficient_;
    double const exponent_;
};

}  // namespace MaterialPropertyLib
