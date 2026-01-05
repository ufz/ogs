// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
struct SecondaryVariable
{
    SecondaryVariable(std::string name_,
                      MeshLib::PropertyVector<double>* value_)
        : name(std::move(name_)), value(value_)
    {
    }

    std::string const name;
    MeshLib::PropertyVector<double>* value;
    static const ItemType item_type = ItemType::SecondaryVariable;
};

struct UserPunch
{
    UserPunch(std::vector<SecondaryVariable>&& secondary_variables_,
              std::vector<std::string>&& statements_)
        : secondary_variables(std::move(secondary_variables_)),
          statements(std::move(statements_))
    {
    }

    void initialize(std::size_t const num_chemical_systems);

    friend std::ostream& operator<<(std::ostream& os,
                                    UserPunch const& user_punch);

    std::vector<SecondaryVariable> secondary_variables;
    std::vector<std::string> const statements;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
