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

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace ChemistryLib
{
class ChemicalSolverInterface
{
public:
    virtual void executeInitialCalculation(
        std::vector<GlobalVector*>& process_solutions) = 0;

    virtual void doWaterChemistryCalculation(
        std::vector<GlobalVector*>& process_solutions, double const dt) = 0;

    virtual ~ChemicalSolverInterface() = default;
};
}  // namespace ChemistryLib
