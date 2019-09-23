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

#include <map>
#include <string>

#include "ThirdParty/iphreeqc/src/src/phreeqcpp/ISolution.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
class Component final : public cxxISolutionComp
{
public:
    explicit Component(std::string const& name) { description = name; }
};

class InitialAqueousSolution final : public cxxISolution
{
public:
    explicit InitialAqueousSolution(
        std::map<std::string, cxxISolutionComp>& components);
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
