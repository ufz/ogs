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
#include <iosfwd>
#include <string>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
struct ReactionRate
{
    ReactionRate(std::string kinetic_reactant_,
                 std::vector<std::string>
                     expression_statements_)
        : kinetic_reactant(std::move(kinetic_reactant_)),
          expression_statements(std::move(expression_statements_))
    {
    }

    friend std::ofstream& operator<<(
        std::ofstream& out, std::vector<ReactionRate> const& reaction_rate);

    std::string const kinetic_reactant;
    std::vector<std::string> const expression_statements;
};

std::vector<ReactionRate> createReactionRates(
    boost::optional<BaseLib::ConfigTree> const& config);
}  // namespace ChemistryLib
