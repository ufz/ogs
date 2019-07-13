/**
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

#include "ChemistryLib/Output.h"
#include "MeshLib/PropertyVector.h"

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
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

    friend std::ostream& operator<<(std::ostream& os,
                                    EquilibriumPhase const& equilibrium_phase);

    std::string const name;
    MeshLib::PropertyVector<double>* amount;
    double const saturation_index;
    static const ItemType item_type = ItemType::EquilibriumPhase;
};
}  // namespace ChemistryLib
