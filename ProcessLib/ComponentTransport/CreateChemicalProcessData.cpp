/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateChemicalProcessData.h"
#include "ChemicalProcessData.h"

#include "ChemistryLib/ChemicalSolverInterface.h"

namespace ProcessLib
{
namespace ComponentTransport
{
std::unique_ptr<ChemicalProcessData> createChemicalProcessData(
    ChemistryLib::ChemicalSolverInterface* const chemical_solver_interface)
{
    if (!chemical_solver_interface)
    {
        return nullptr;
    }

    return std::make_unique<ChemicalProcessData>(
        chemical_solver_interface->chemical_system_index_map);
}
}  // namespace ComponentTransport
}  // namespace ProcessLib
