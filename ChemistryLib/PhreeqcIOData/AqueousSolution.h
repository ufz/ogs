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

#include "Output.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct Component
{
    explicit Component(std::string name_) : name(std::move(name_)) {}

    std::string const name;
    double amount = std::numeric_limits<double>::quiet_NaN();
    static const ItemType item_type = ItemType::Component;
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

    friend std::ostream& operator<<(std::ostream& os,
                                    AqueousSolution const& aqueous_solution);

    double temperature;
    double pressure;
    double pH = std::numeric_limits<double>::quiet_NaN();
    double pe;
    std::vector<Component> components;
    MeansOfAdjustingCharge const means_of_adjusting_charge;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
