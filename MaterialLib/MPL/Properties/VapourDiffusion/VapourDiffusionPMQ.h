/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
 * The model was presented in \cite chau2005simulation
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
 *  In the PMQ type, \f$D_{vr}\f$ takes the form of \cite chau2005simulation
 *   \f[
 *     D_{vr}=0.66 \phi (1 - S_L )^2,
 *   \f]
 *    with \f$\phi\f$, the porosity, \f$S_L\f$, the water saturation, \f$0.66\f$
 *     the tortuosity constant.
 *
 */

class VapourDiffusionPMQ final : public Property
{
public:
    explicit VapourDiffusionPMQ(std::string name) { name_ = std::move(name); }

    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'VapourDiffusionPMQ' is "
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
    static double constexpr tortuosity_ = 0.66;
};

}  // namespace MaterialPropertyLib
