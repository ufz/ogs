/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
                        MeshLib::PropertyVector<double>* molality_,
                        MeshLib::PropertyVector<double>* volume_fraction_,
                        MeshLib::PropertyVector<double>* mesh_prop_molality_,
                        double saturation_index_)
        : name(std::move(name_)),
          molality(molality_),
          volume_fraction(volume_fraction_),
          mesh_prop_molality(mesh_prop_molality_),
          saturation_index(saturation_index_)
    {
    }

    void print(std::ostream& os, std::size_t const chemical_system_id) const;

    std::string const name;
    MeshLib::PropertyVector<double>* molality;
    MeshLib::PropertyVector<double>* volume_fraction;
    MeshLib::PropertyVector<double>* mesh_prop_molality;
    double const saturation_index;
    static const ItemType item_type = ItemType::EquilibriumReactant;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
