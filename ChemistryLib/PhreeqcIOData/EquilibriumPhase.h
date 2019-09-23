/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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
struct EquilibriumPhase
{
    EquilibriumPhase(std::string name_,
                     MeshLib::PropertyVector<double>* amount_,
                     double saturation_index_)
        : name(std::move(name_)),
          amount(amount_),
          saturation_index(saturation_index_)
    {
    }

    void print(std::ostream& os, std::size_t const chemical_system_id) const;

    std::string const name;
    MeshLib::PropertyVector<double>* amount;
    double const saturation_index;
    static const ItemType item_type = ItemType::EquilibriumPhase;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
