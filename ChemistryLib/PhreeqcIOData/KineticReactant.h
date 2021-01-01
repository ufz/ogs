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

#include <boost/optional/optional.hpp>
#include <iosfwd>
#include <string>
#include <vector>

#include "MeshLib/PropertyVector.h"
#include "Output.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct KineticReactant
{
    KineticReactant(std::string name_,
                    std::string chemical_formula_,
                    MeshLib::PropertyVector<double>* molality_,
                    MeshLib::PropertyVector<double>* mesh_prop_molality_,
                    std::vector<double>&& parameters_,
                    bool const fix_amount_)
        : name(std::move(name_)),
          chemical_formula(std::move(chemical_formula_)),
          molality(molality_),
          mesh_prop_molality(mesh_prop_molality_),
          parameters(std::move(parameters_)),
          fix_amount(fix_amount_)
    {
    }

    void print(std::ostream& os, std::size_t const global_id) const;

    std::string const name;
    std::string const chemical_formula;
    MeshLib::PropertyVector<double>* molality;
    MeshLib::PropertyVector<double>* mesh_prop_molality;
    std::vector<double> const parameters;
    bool const fix_amount;
    static const ItemType item_type = ItemType::KineticReactant;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
