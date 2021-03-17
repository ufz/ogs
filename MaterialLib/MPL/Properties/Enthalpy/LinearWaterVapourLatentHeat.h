/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 16, 2021, 10:03 AM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Phase;

/**
 * \brief An empirical function for the latent heat of vaporization
 *      of liquid water, which is given by \cite saito2006numerical p.786f
 *  \f[
 *     L_w(T)=2.501 \cdot 10^6 - 2369.2  (T - 273.15),\,\text{[J/kg]}.
 *  \f]
 *
 *  A quite simular formula is present on page 81 of \cite bittelli2015soil.
 *
 *  The linear expressions of the latent heat of vaporization of liquid water
 *  can be found in some very early references, which are mentioned in
 *  \cite harrison1965fundamental (page 79).
 *
 *  The function is used in the energy balance equation for partially saturated
 *  zone in porous media \cite de1958simultaneous.
 */
class LinearWaterVapourLatentHeat final : public Property
{
public:
    explicit LinearWaterVapourLatentHeat(std::string name)
    {
        name_ = std::move(name);
    }

    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'LinearWaterVapourLatentHeat' is "
                "implemented on the 'phase' scale only.");
        }
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};

}  // namespace MaterialPropertyLib
