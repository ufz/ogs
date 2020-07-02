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

#include <memory>

namespace ChemistryLib
{
class ChemicalSolverInterface;
}

namespace ProcessLib
{
namespace ComponentTransport
{
struct ChemicalProcessData;

std::unique_ptr<ChemicalProcessData> createChemicalProcessData(
    ChemistryLib::ChemicalSolverInterface* chemical_solver_interface);
}  // namespace ComponentTransport
}  // namespace ProcessLib
