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

#include "ChemistryLib/Output.h"

namespace ChemistryLib
{
struct KineticReactant
{
    KineticReactant(std::string name_,
                    double amount_,
                    boost::optional<std::vector<double>>
                        parameters_)
        : name(std::move(name_)),
          amount(amount_),
          parameters(parameters_ ? std::move(*parameters_)
                                 : std::vector<double>{})
    {
    }

    friend std::ofstream& operator<<(
        std::ofstream& out,
        std::vector<KineticReactant> const& kinetic_reactants);

    std::string const name;
    double amount;
    std::vector<double> const parameters;
    static const ItemType item_type = ItemType::KineticReactant;
};
}  // namespace ChemistryLib
