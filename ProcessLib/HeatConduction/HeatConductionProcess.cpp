// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "HeatConductionProcess.h"

#include <cassert>

#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#include "MathLib/LinAlg/FinalizeVectorAssembly.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Utils/ComputeResiduum.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace HeatConduction
{
HeatConductionProcess::HeatConductionProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    HeatConductionProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const is_linear,
    bool const ls_compute_only_upon_timestep_change)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      AssemblyMixin<HeatConductionProcess>{*_jacobian_assembler, is_linear,
                                           true /*use_monolithic_scheme*/},
      _process_data(std::move(process_data)),
      _ls_compute_only_upon_timestep_change{
          ls_compute_only_upon_timestep_change}
{
    // For numerical Jacobian assembler
    this->_jacobian_assembler->setNonDeformationComponentIDs(
        {0} /* only one variable: temperature */);

    if (ls_compute_only_upon_timestep_change)
    {
        // TODO move this feature to some common location for all processes.
        if (!is_linear)
        {
            OGS_FATAL(
                "Using the linear solver compute() method only upon timestep "
                "change only makes sense for linear model equations.");
        }

        WARN(
            "You specified that the HeatConduction linear solver will do "
            "the compute() step only upon timestep change. This is an expert "
            "option. It is your responsibility to ensure that "
            "the conditions for the correct use of this feature are met! "
            "Otherwise OGS might compute garbage without being recognized. "
            "There is no safety net!");

        WARN(
            "You specified that the HeatConduction linear solver will do "
            "the compute() step only upon timestep change. This option will "
            "only be used by the Picard non-linear solver. The Newton-Raphson "
            "non-linear solver will silently ignore this setting, i.e., it "
            "won't do any harm, there, but you won't observe the speedup you "
            "probably expect.");
    }
}

void HeatConductionProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    int const mesh_space_dimension = _process_data.mesh_space_dimension;

    ProcessLib::createLocalAssemblers<LocalAssemblerData>(
        mesh_space_dimension, mesh.getElements(), dof_table, local_assemblers_,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        _process_data);

    _secondary_variables.addSecondaryVariable(
        "heat_flux",
        makeExtrapolator(
            mesh_space_dimension, getExtrapolator(), local_assemblers_,
            &HeatConductionLocalAssemblerInterface::getIntPtHeatFlux));
}

std::vector<std::vector<std::string>>
HeatConductionProcess::initializeAssemblyOnSubmeshes(
    std::vector<std::reference_wrapper<MeshLib::Mesh>> const& meshes)
{
    DBUG("HeatConductionProcess initializeSubmeshOutput().");

    std::vector<std::vector<std::string>> residuum_names{{"HeatFlowRate"}};

    AssemblyMixin<HeatConductionProcess>::initializeAssemblyOnSubmeshes(
        meshes, residuum_names);

    return residuum_names;
}

void HeatConductionProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble HeatConductionProcess.");

    AssemblyMixin<HeatConductionProcess>::assemble(t, dt, x, x_prev, process_id,
                                                   M, K, b);
}

void HeatConductionProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian HeatConductionProcess.");

    AssemblyMixin<HeatConductionProcess>::assembleWithJacobian(
        t, dt, x, x_prev, process_id, b, Jac);
}

void HeatConductionProcess::computeSecondaryVariableConcrete(
    double const t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& x_prev, int const process_id)
{
    DBUG("Compute heat flux for HeatConductionProcess.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    std::generate_n(std::back_inserter(dof_tables), x.size(),
                    [&]() { return _local_to_global_index_map.get(); });

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &HeatConductionLocalAssemblerInterface::computeSecondaryVariable,
        local_assemblers_, getActiveElementIDs(), dof_tables, t, dt, x, x_prev,
        process_id);
}

void HeatConductionProcess::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& /*x*/,
    const double /*t*/,
    const double /*dt*/,
    const int /*process_id*/)
{
    AssemblyMixin<HeatConductionProcess>::updateActiveElements();
}

void HeatConductionProcess::preOutputConcreteProcess(
    const double t,
    double const dt,
    std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev,
    int const process_id)
{
    AssemblyMixin<HeatConductionProcess>::preOutput(t, dt, x, x_prev,
                                                    process_id);
}
}  // namespace HeatConduction
}  // namespace ProcessLib
