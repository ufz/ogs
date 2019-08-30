/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ThirdParty/iphreeqc/src/src/phreeqcpp/cxxKinetics.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
class KineticReactant final : public cxxKineticsComp
{
public:
    KineticReactant(std::string const& name, double const initial_amount);
};

class Kinetics final : public cxxKinetics
{
public:
    explicit Kinetics(std::vector<KineticReactant> const& kinetic_reactants);
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
