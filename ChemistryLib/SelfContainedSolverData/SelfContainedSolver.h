/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Sparse>
#include <memory>

#include "ChemicalReaction.h"
#include "ChemistryLib/ChemicalSolverInterface.h"

namespace ChemistryLib
{
namespace SelfContainedSolverData
{
class SelfContainedSolver final : public ChemicalSolverInterface
{
public:
    SelfContainedSolver(MeshLib::Mesh const& mesh,
                        GlobalLinearSolver& linear_solver,
                        Eigen::SparseMatrix<double>
                            stoichiometric_matrix,
                        std::vector<std::unique_ptr<ChemicalReaction>>
                            chemical_reactions)
        : ChemicalSolverInterface(mesh, linear_solver),
          _stoichiometric_matrix(stoichiometric_matrix),
          _chemical_reactions(std::move(chemical_reactions))
    {
    }

    Eigen::SparseMatrix<double> const* getStoichiometricMatrix() const override
    {
        return &_stoichiometric_matrix;
    }

    double getKineticPrefactor(std::size_t reaction_id) const override
    {
        return _chemical_reactions[reaction_id]->getKineticPrefactor();
    }

private:
    Eigen::SparseMatrix<double> _stoichiometric_matrix;
    std::vector<std::unique_ptr<ChemicalReaction>> _chemical_reactions;
};
}  // namespace SelfContainedSolverData
}  // namespace ChemistryLib
