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

#include <phreeqcpp/cxxKinetics.h>

#include <string>
#include <vector>

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
class KineticReactant final : private cxxKineticsComp
{
public:
    KineticReactant(std::string const& name, double const initial_amount);

    cxxKineticsComp const* castToBaseClass() const
    {
        return static_cast<cxxKineticsComp const*>(this);
    }
};

class Kinetics final : private cxxKinetics
{
public:
    explicit Kinetics(std::vector<KineticReactant> const& kinetic_reactants);

    void setChemicalSystemID(std::size_t const chemical_system_id)
    {
        Set_n_user_both(chemical_system_id);
    }

    cxxKinetics const* castToBaseClass() const
    {
        return static_cast<cxxKinetics const*>(this);
    }
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
