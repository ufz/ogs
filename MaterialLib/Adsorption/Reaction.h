/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>

namespace BaseLib { class ConfigTree; }

namespace Adsorption
{

class Reaction
{
public:
    static std::unique_ptr<Reaction> newInstance(
        BaseLib::ConfigTree const& conf);

    virtual double getEnthalpy(const double p_Ads, const double T_Ads, const double M_Ads) const = 0;
    virtual double getReactionRate(const double p_Ads, const double T_Ads,
                                     const double M_Ads, const double loading) const = 0;

    // TODO get rid of
    virtual double getEquilibriumLoading(const double /*unused*/,
                                         const double /*unused*/,
                                         const double /*unused*/) const
    {
        return -1.0;
    }

    virtual ~Reaction() = default;
};

}  // namespace Adsorption
