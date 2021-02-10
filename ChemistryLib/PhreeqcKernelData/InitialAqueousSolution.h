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

#include <phreeqcpp/ISolution.h>

#include <map>
#include <string>

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
