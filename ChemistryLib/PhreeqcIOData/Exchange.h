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

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct Exchange
{
    Exchange(std::string name_, double const molality_)
        : name(std::move(name_)), molality(molality_)
    {
    }

    friend std::ostream& operator<<(std::ostream& os,
                                    Exchange const& exchange_assemblage);

    std::string const name;
    double const molality;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
