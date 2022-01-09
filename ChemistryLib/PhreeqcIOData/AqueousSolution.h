/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "Output.h"

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ChemistryLib
{
enum class ChargeBalance;

namespace PhreeqcIOData
{
struct Component
{
    explicit Component(std::string name_, std::string chemical_formula_)
        : name(std::move(name_)), chemical_formula(std::move(chemical_formula_))
    {
    }

    std::string const name;
    std::string const chemical_formula;
    std::unique_ptr<GlobalVector> amount;
    static const ItemType item_type = ItemType::Component;
};

struct AqueousSolution
{
    AqueousSolution(double temperature_, double pressure_,
                    MeshLib::PropertyVector<double>* pe_, double const pe0_,
                    std::vector<Component>&& components_,
                    ChargeBalance charge_balance_)
        : temperature(temperature_),
          pressure(pressure_),
          pe(pe_),
          pe0(pe0_),
          components(std::move(components_)),
          charge_balance(charge_balance_)
    {
    }

    void print(std::ostream& os, std::size_t const chemical_system_id) const;

    double const temperature;
    double const pressure;
    std::unique_ptr<GlobalVector> pH;
    MeshLib::PropertyVector<double>* pe;
    double const pe0;
    std::vector<Component> components;
    ChargeBalance const charge_balance;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
