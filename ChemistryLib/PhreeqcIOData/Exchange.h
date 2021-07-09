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
