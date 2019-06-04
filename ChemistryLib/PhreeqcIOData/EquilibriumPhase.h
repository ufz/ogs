/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/optional/optional.hpp>
#include <string>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
struct EquilibriumPhase
{
    EquilibriumPhase(std::string name_,
                     double amount_,
                     double saturation_index_)
        : name(std::move(name_)),
          amount(amount_),
          saturation_index(saturation_index_)
    {
    }

    std::string const name;
    double amount;
    double const saturation_index;
};

std::vector<EquilibriumPhase> createEquilibriumPhases(
    boost::optional<BaseLib::ConfigTree> const& config);
}  // namespace ChemistryLib
