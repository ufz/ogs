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

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace MaterialPropertyLib
{
class Medium;
}

namespace ParameterLib
{
class SpatialPosition;
}

namespace ChemistryLib
{
class ChemicalSolverInterface
{
public:
    virtual void initialize() {}

    virtual void initializeChemicalSystemConcrete(
        GlobalIndexType const& /*chemical_system_id*/,
        MaterialPropertyLib::Medium const* /*medium*/,
        ParameterLib::SpatialPosition const& /*pos*/,
        double const /*t*/)
    {
    }

    virtual void executeInitialCalculation() = 0;

    virtual void doWaterChemistryCalculation(double const dt) = 0;

    virtual std::vector<GlobalVector*> getIntPtProcessSolutions() const = 0;

    virtual std::vector<std::string> const getComponentList() const
    {
        return {};
    }

    virtual void computeSecondaryVariable(
        std::size_t const /*ele_id*/,
        std::vector<GlobalIndexType> const& /*chemical_system_indices*/)
    {
    }

    virtual ~ChemicalSolverInterface() = default;

public:
    std::vector<GlobalIndexType> chemical_system_index_map;
};
}  // namespace ChemistryLib
