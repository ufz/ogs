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
#include "MeshLib/PropertyVector.h"

namespace ChemistryLib
{
struct KineticReactant
{
    KineticReactant(std::string&& name_,
                    std::string&& chemical_formula_,
                    MeshLib::PropertyVector<double>* amount_,
                    std::vector<double>&& parameters_)
        : name(std::move(name_)),
          chemical_formula(std::move(chemical_formula_)),
          amount(amount_),
          parameters(std::move(parameters_))
    {
    }

    void print(std::ostream& os, std::size_t const chemical_system_id) const;

    std::string const name;
    std::string const chemical_formula;
    MeshLib::PropertyVector<double>* amount;
    std::vector<double> const parameters;
    static const ItemType item_type = ItemType::KineticReactant;
};
}  // namespace ChemistryLib
