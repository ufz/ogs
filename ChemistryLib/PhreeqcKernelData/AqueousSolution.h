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

// clang-format off
#include <phreeqcpp/common/phrqtype.h>
#include <phreeqcpp/Solution.h>
// clang-format on

#include <memory>

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
class AqueousSolution final : private cxxSolution
{
public:
    AqueousSolution(double const temperature, double const pressure,
                    double const pe_value,
                    cxxISolution const& initial_aqueous_solution);

    std::unique_ptr<cxxISolution const> getInitialAqueousSolution() const
    {
        return std::make_unique<cxxISolution const>(*Get_initial_data());
    }

    void setChemicalSystemID(std::size_t const chemical_system_id)
    {
        Set_n_user_both(chemical_system_id);
    }

    cxxSolution const* castToBaseClass() const
    {
        return static_cast<cxxSolution const*>(this);
    }

    std::unique_ptr<cxxSolution const> castToBaseClassNoninitialized()
    {
        Destroy_initial_data();
        return std::make_unique<cxxSolution const>(*castToBaseClass());
    }
};
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
