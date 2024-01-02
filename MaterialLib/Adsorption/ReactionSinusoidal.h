/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Reaction.h"

namespace Adsorption
{

class ReactionSinusoidal final : public Reaction
{
public:
    explicit ReactionSinusoidal(BaseLib::ConfigTree const& conf);

    double getEnthalpy(const double /*p_Ads*/, const double /*T_Ads*/,
                       const double /*M_Ads*/) const override
    {
        return _enthalpy;
    }

    double getReactionRate(const double /*p_Ads*/, const double /*T_Ads*/,
                           const double /*M_Ads*/,
                           const double /*loading*/) const override;

private:
    double _enthalpy;
};

}  // namespace Adsorption
