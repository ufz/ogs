/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Sparse>

#include "MaterialLib/MPL/VariableType.h"
#include "MathLib/LinAlg/GlobalLinearSolverType.h"
#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MeshLib/Mesh.h"

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
    ChemicalSolverInterface(MeshLib::Mesh const& mesh,
                            GlobalLinearSolver& linear_solver_)
        : _mesh(mesh), linear_solver(linear_solver_)
    {
        auto const* const bulk_element_ids = bulkElementIDs(_mesh);
        if (bulk_element_ids == nullptr)
        {
            OGS_FATAL(
                "The 'bulk_element_ids' property does not exist on the mesh "
                "{:s}.",
                _mesh.getName());
        }
        active_element_ids_.assign(bulk_element_ids->begin(),
                                   bulk_element_ids->end());
    }

    std::vector<std::size_t> const& activeElementIDs() const
    {
        return active_element_ids_;
    }

    virtual void initialize() {}

    virtual void initializeChemicalSystemConcrete(
        std::vector<double> const& /*concentrations*/,
        GlobalIndexType const& /*chemical_system_id*/,
        MaterialPropertyLib::Medium const& /*medium*/,
        ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/)
    {
    }

    virtual void setChemicalSystemConcrete(
        std::vector<double> const& /*concentrations*/,
        GlobalIndexType const& /*chemical_system_id*/,
        MaterialPropertyLib::Medium const* /*medium*/,
        MaterialPropertyLib::VariableArray const& /*vars*/,
        ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
        double const /*dt*/)
    {
    }

    virtual void setAqueousSolutionsPrevFromDumpFile() {}

    virtual void executeSpeciationCalculation([[maybe_unused]] double const dt)
    {
    }

    virtual double getConcentration(
        int const /*component_id*/,
        GlobalIndexType const /*chemical_system_id*/) const
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    virtual std::vector<std::string> const getComponentList() const
    {
        return {};
    }

    virtual Eigen::SparseMatrix<double> const* getStoichiometricMatrix() const
    {
        return nullptr;
    }

    virtual double getKineticPrefactor(
        [[maybe_unused]] std::size_t reaction_id) const
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    virtual void updateVolumeFractionPostReaction(
        GlobalIndexType const& /*chemical_system_id*/,
        MaterialPropertyLib::Medium const& /*medium*/,
        ParameterLib::SpatialPosition const& /*pos*/, double const /*porosity*/,
        double const /*t*/, double const /*dt*/)
    {
    }

    virtual void updatePorosityPostReaction(
        GlobalIndexType const& /*chemical_system_id*/,
        MaterialPropertyLib::Medium const& /*medium*/,
        double& /*porosity*/)
    {
    }

    virtual void computeSecondaryVariable(
        std::size_t const /*ele_id*/,
        std::vector<GlobalIndexType> const& /*chemical_system_indices*/)
    {
    }

    virtual ~ChemicalSolverInterface() = default;

public:
    std::vector<GlobalIndexType> chemical_system_index_map;
    MeshLib::Mesh const& _mesh;
    /// specify the linear solver used to solve the linearized reaction
    /// equation.
    GlobalLinearSolver& linear_solver;

private:
    std::vector<std::size_t> active_element_ids_;
};
}  // namespace ChemistryLib
