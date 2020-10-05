/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "MeshLib/PropertyVector.h"
#include "Output.h"

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct EquilibriumReactant
{
    EquilibriumReactant(std::string name_,
                        MeshLib::PropertyVector<double>* amount_,
                        MeshLib::PropertyVector<double>* amount_avg_,
                        double const initial_amount_,
                        double saturation_index_)
        : name(std::move(name_)),
          amount(amount_),
          amount_avg(amount_avg_),
          initial_amount(initial_amount_),
          saturation_index(saturation_index_)
    {
    }

    void print(std::ostream& os, std::size_t const global_id) const;

    std::string const name;
    MeshLib::PropertyVector<double>* amount;
    MeshLib::PropertyVector<double>* amount_avg;
    double const initial_amount;
    double const saturation_index;
    static const ItemType item_type = ItemType::EquilibriumReactant;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
