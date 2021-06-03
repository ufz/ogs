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
struct ExchangeSite
{
    ExchangeSite(std::string ion_exchanging_species_, double const molality_)
        : ion_exchanging_species(std::move(ion_exchanging_species_)),
          molality(molality_)
    {
    }

    std::string const ion_exchanging_species;
    double const molality;
};

std::ostream& operator<<(std::ostream& os, ExchangeSite const& exchange_site);

}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
