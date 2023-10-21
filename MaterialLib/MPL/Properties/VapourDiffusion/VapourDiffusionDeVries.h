/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 30, 2023, 1:49 PM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Phase;

/**
 * \brief DeVries type Vapour diffusion
 *
 *  The model was presented in \cite DeVries1957.
 *
 *  The vapour diffusion can be described by
 *  \f[
 *     D_v=5.9\cdot 10^{-6} \left(\frac{T, \text{in K}}{p, \text{in
 * Pa}}\right)^{2.3} D_{vr}, \f] where \f$D_{vr}\f$ is the the relative
 * diffusion coefficient, and \f$T\f$ is the temperature.
 *
 */
class VapourDiffusionDeVries final : public Property
{
public:
    VapourDiffusionDeVries(std::string name,
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
                "The property 'VapourDiffusionDeVries' is "
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
