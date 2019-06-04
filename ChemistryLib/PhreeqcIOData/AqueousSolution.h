/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "BaseLib/Error.h"

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
struct Component
{
    Component(std::string name_) : name(std::move(name_)) {}

    std::string const name;
    double amount;
};

enum class MeansOfAdjustingCharge
{
    pH,
    pe,
    Unspecified
};

struct AqueousSolution
{
    AqueousSolution(double temperature_,
                    double pressure_,
                    double pe_,
                    std::vector<Component>&& components_,
                    MeansOfAdjustingCharge means_of_adjusting_charge_)
        : temperature(temperature_),
          pressure(pressure_),
          pe(pe_),
          components(std::move(components_)),
          means_of_adjusting_charge(means_of_adjusting_charge_)
    {
    }

    friend std::ofstream& operator<<(std::ofstream& out,
                                     AqueousSolution const& aqueous_solution);

    double temperature;
    double pressure;
    double pH;
    double pe;
    std::vector<Component> components;
    MeansOfAdjustingCharge const means_of_adjusting_charge;
};

AqueousSolution createAqueousSolution(BaseLib::ConfigTree const& config);
}  // namespace ChemistryLib
