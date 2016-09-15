/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_SMALLDEFORMATIONPROCESS_H_
#define PROCESS_LIB_SMALLDEFORMATIONPROCESS_H_

#include <cassert>

#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

#include "SmallDeformationFEM.h"
#include "SmallDeformationProcessData.h"

namespace ProcessLib
{
namespace SmallDeformation
{
template <int DisplacementDim>
class SmallDeformationProcess final : public Process
{
    using Base = Process;

public:
    SmallDeformationProcess(
        MeshLib::Mesh& mesh,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        SmallDeformationProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller)
        : Process(mesh, parameters, std::move(process_variables),
                  std::move(secondary_variables),
                  std::move(named_function_caller)),
          _process_data(std::move(process_data))
    {
    }

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

private:
    using LocalAssemblerInterface = SmallDeformationLocalAssemblerInterface;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override
    {
        ProcessLib::SmallDeformation::createLocalAssemblers<DisplacementDim,
                                                            LocalAssemblerData>(
            mesh.getDimension(), mesh.getElements(), dof_table,
            integration_order, _local_assemblers, _process_data);

        // TODO move the two data members somewhere else.
        // for extrapolation of secondary variables
        std::vector<std::unique_ptr<MeshLib::MeshSubsets>>
            all_mesh_subsets_single_component;
        all_mesh_subsets_single_component.emplace_back(
            new MeshLib::MeshSubsets(_mesh_subset_all_nodes.get()));
        _local_to_global_index_map_single_component.reset(
            new NumLib::LocalToGlobalIndexMap(
                std::move(all_mesh_subsets_single_component),
                // by location order is needed for output
                NumLib::ComponentOrder::BY_LOCATION));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xx", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXX));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_yy", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtSigmaYY));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xy", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXY));
    }

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override
    {
        DBUG("Assemble SmallDeformationProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberOnDereferenced(
            &SmallDeformationLocalAssemblerInterface::assemble,
            _local_assemblers, *_local_to_global_index_map, t, x, M, K, b);
    }

    void assembleJacobianConcreteProcess(const double t, GlobalVector const& x,
                                         GlobalVector const& /*xdot*/,
                                         const double /*dxdot_dx*/,
                                         GlobalMatrix const& /*M*/,
                                         const double /*dx_dx*/,
                                         GlobalMatrix const& /*K*/,
                                         GlobalMatrix& Jac) override
    {
        DBUG("AssembleJacobian SmallDeformationProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberOnDereferenced(
            &SmallDeformationLocalAssemblerInterface::assembleJacobian,
            _local_assemblers, *_local_to_global_index_map, t, x, Jac);
    }

    void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                                    double const dt) override
    {
        DBUG("PreTimestep SmallDeformationProcess.");

        _process_data.dt = dt;
        _process_data.t = t;

        GlobalExecutor::executeMemberOnDereferenced(
            &SmallDeformationLocalAssemblerInterface::preTimestep,
            _local_assemblers, *_local_to_global_index_map, x, t, dt);
    }

    void postTimestepConcreteProcess(GlobalVector const& x) override
    {
        DBUG("PostTimestep SmallDeformationProcess.");

        GlobalExecutor::executeMemberOnDereferenced(
            &SmallDeformationLocalAssemblerInterface::postTimestep,
            _local_assemblers, *_local_to_global_index_map, x);
    }

private:
    SmallDeformationProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;
};

}  // namespace SmallDeformation
}  // namespace ProcessLib

#endif  // PROCESS_LIB_SMALLDEFORMATIONPROCESS_H_
