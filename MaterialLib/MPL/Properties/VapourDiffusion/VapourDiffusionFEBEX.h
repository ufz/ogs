/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
 * 	D_v=2.16\cdot 10^{-5} \left(\frac{T}{273.15}\right)^{1.8}
 * 	D_{vr},
 *  \f]
 *  where \f$D_{vr}\f$ is the the relative diffusion coefficient,
 *  and \f$T\f$ is the temperature.
 *
 *  In the FEBEX type, \f$D_{vr}\f$ takes the form of \cite Rutquist2007TaskA1
 *   \f[
 *     D_{vr}=\tau \phi (1 - S ),
 *   \f]
 *    with \f$\phi\f$, the porosity, \f$S\f$, the water saturation,
 *    and  \f$\tau\f$ the tortuosity.
 *
 */
class VapourDiffusionFEBEX final : public Property
{
public:
    VapourDiffusionFEBEX(std::string name, double const tortuosity)
        : tortuosity_(tortuosity)
    {
        name_ = std::move(name);
    }

    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'VapourDiffusionFEBEX' is "
                "implemented on the 'phase' scale only.");
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
    double const tortuosity_;
};

}  // namespace MaterialPropertyLib
