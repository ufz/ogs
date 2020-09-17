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

#include "BoundaryCondition.h"
#include "ParameterLib/MeshNodeParameter.h"

namespace ProcessLib
{
/// The SolutionDependentDirichletBoundaryCondition belongs to the
/// Dirichlet-type boundary condition.
///
/// This class is a special category of Dirichlet boundary condition,
/// applied in the situation where the value assigned for the boundary condition
/// is dependent on the process solution of the last time step. This particular
/// boundary condition is widely used in the reactive transport problems and has
/// the potential to be used in other processes.
class SolutionDependentDirichletBoundaryCondition final
    : public BoundaryCondition
{
public:
    SolutionDependentDirichletBoundaryCondition(
        ParameterLib::Parameter<double> const& parameter,
        MeshLib::Mesh const& bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id);

    void getEssentialBCValues(
        double const t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

    /// Renchao: The original idea to place the update of the boundary condition
    /// value at the preTimestep stage. The update could be achieved within the
    /// class function "Process::preTimestep". Whereas, I find it not doable to
    /// implement in this way. The class function "Process::preTimestep" is
    /// called in a row by the functions "preTimestepForAllProcesses" and
    /// "TimeLoop::outputSolutions". These two functions are called when
    /// initializing and subsequently looping over the "TimeLoop". Actually, it
    /// is not intended to make the boundary condition value updated in the
    /// first loop. Instead, the update is intended to start from the second
    /// loop. For these two reasons, I think it more proper to do the
    /// implementation at the postTimestep stage.
    void postTimestep(double const /*t*/,
                      std::vector<GlobalVector*> const& x,
                      int const process_id) override;

private:
    MeshLib::Mesh const& _bc_mesh;
    int const _variable_id;
    int const _component_id;
    std::unique_ptr<NumLib::LocalToGlobalIndexMap const> _dof_table_boundary;
    std::unique_ptr<ParameterLib::MeshNodeParameter<double>> _parameter;
};

std::unique_ptr<SolutionDependentDirichletBoundaryCondition>
createSolutionDependentDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace ProcessLib
