// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include "MeshLib/PropertyVector.h"

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct ExchangeSite
{
    ExchangeSite(std::string name_, MeshLib::PropertyVector<double>* molality_)
        : name(std::move(name_)), molality(molality_)
    {
    }

    void print(std::ostream& os, std::size_t const chemical_system_id) const;

    std::string const name;
    MeshLib::PropertyVector<double>* molality;
};

}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
