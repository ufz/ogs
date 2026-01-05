// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
