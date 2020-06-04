/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "Output.h"

namespace ChemistryLib
{
enum class ChargeBalance;

namespace PhreeqcIOData
{
struct Component
{
    explicit Component(std::string name_,
                       std::size_t const num_chemical_systems_)
        : name(std::move(name_)),
          amount(MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
              num_chemical_systems_))
    {
    }

    std::string const name;
    std::unique_ptr<GlobalVector> amount;
    static const ItemType item_type = ItemType::Component;
};

struct AqueousSolution
{
    AqueousSolution(double temperature_,
                    double pressure_,
                    double pe_,
                    std::vector<Component>&& components_,
                    ChargeBalance charge_balance_,
                    std::size_t const num_chemical_systems_)
        : temperature(temperature_),
          pressure(pressure_),
          pH(MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
              num_chemical_systems_)),
          pe(pe_),
          components(std::move(components_)),
          charge_balance(charge_balance_)
    {
    }

    friend std::ostream& operator<<(std::ostream& os,
                                    AqueousSolution const& aqueous_solution);

    double temperature;
    double pressure;
    double pe;
    std::unique_ptr<GlobalVector> pH;
    std::vector<Component> components;
    ChargeBalance const charge_balance;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
