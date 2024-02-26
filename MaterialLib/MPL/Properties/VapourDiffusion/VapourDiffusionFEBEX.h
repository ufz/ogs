/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 5, 2021, 3:49 PM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Phase;

/**
 * \brief FEBEX type Vapour diffusion
 *
 *  The model was presented in \cite Rutquist2007TaskA1.
 *
 *  The vapour diffusion can be described by
 *   \cite moldrup1997modeling, \cite moldrup2000predicting,
 *  \cite chau2005simulation
 *  \f[
 *     D_v=D_0 \left(\frac{T}{273.15}\right)^{n}
 *     D_{vr},
 *  \f]
 *  where \f$D_{0}\f$ is the base diffusion coefficient with default value
 *  \f$2.16\cdot 10^{-5}\f$ \f${\text m}^2
 * \text{Pa}/(\text{s}\text{K}^{n})\f$,
 *  \f$n\f$ is the exponent with default value 1.8,
 *  \f$D_{vr}\f$ is the the relative
 * diffusion coefficient, and \f$T\f$ is the temperature.
 *
 *  In the FEBEX type, \f$D_{vr}\f$ takes the form of \cite Rutquist2007TaskA1
 *   \f[
 *     D_{vr}=\tau \phi (1-S_L),
 *   \f]
 *    with \f$\phi\f$, the porosity, \f$S_L\f$, the liquid saturation,
 *    and  \f$\tau\f$, the tortuosity.
 *
 *  Note: In order to maintain consistency with the implementation of the
 *  computations of other vapor-related parameters, \f$ \phi (1 - S_L )\f$ is
 *  removed from the implementation for this class and is multiplied back in
 *  the local assembler.
 */
class VapourDiffusionFEBEX final : public Property
{
public:
    VapourDiffusionFEBEX(std::string name,
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
                "The property 'VapourDiffusionFEBEX' is "
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
