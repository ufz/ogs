/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

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
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        SmallDeformationProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller)
        : Process(mesh, std::move(jacobian_assembler), parameters,
                  integration_order, std::move(process_variables),
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
            _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
            _process_data);

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
            "sigma_zz", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtSigmaZZ));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xy", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXY));

        if (DisplacementDim == 3) {
            Base::_secondary_variables.addSecondaryVariable(
                "sigma_xz", 1,
                makeExtrapolator(
                    getExtrapolator(), _local_assemblers,
                    &SmallDeformationLocalAssemblerInterface::getIntPtSigmaXZ));

            Base::_secondary_variables.addSecondaryVariable(
                "sigma_yz", 1,
                makeExtrapolator(
                    getExtrapolator(), _local_assemblers,
                    &SmallDeformationLocalAssemblerInterface::getIntPtSigmaYZ));
        }

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xx", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonXX));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_yy", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonYY));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_zz", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonZZ));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xy", 1,
            makeExtrapolator(
                getExtrapolator(), _local_assemblers,
                &SmallDeformationLocalAssemblerInterface::getIntPtEpsilonXY));

        Base::integration_point_writer = [this](
            MeshLib::PropertyVector<char>& output,
            MeshLib::PropertyVector<std::size_t>& offsets) {
            return writeIntegrationPointData(output, offsets);
        };
    }

    std::size_t writeIntegrationPointData(MeshLib::PropertyVector<char>& output,
            MeshLib::PropertyVector<std::size_t>& offsets)
    {
        output.clear();
        offsets.clear();
        std::vector<char> local_data;
        std::size_t offset = 0;
        for (auto& la : _local_assemblers)
        {
            offsets.push_back(offset);
            std::size_t const local_offset =
                la->writeIntegrationPointData(local_data);
            std::copy_n(std::begin(local_data), local_offset,
                        std::back_inserter(output));
            offset += local_offset;
        }
        return offset;
    }

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b,
                                 StaggeredCouplingTerm const&
                                 coupling_term) override
    {
        DBUG("Assemble SmallDeformationProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assemble,
            _local_assemblers, *_local_to_global_index_map, t, x, M, K, b,
            coupling_term);
    }

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac,
        StaggeredCouplingTerm const& coupling_term) override
    {
        DBUG("AssembleWithJacobian SmallDeformationProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
            _local_assemblers, *_local_to_global_index_map, t, x, xdot,
            dxdot_dx, dx_dx, M, K, b, Jac, coupling_term);
    }

    void setInitialConditionsConcreteProcess(double const t,
                                             GlobalVector const& x) override
    {
        DBUG("SetInitialConditions SmallDeformationProcess.");

        if (!_mesh.getProperties().hasPropertyVector("integration_point_data"))
            return;
        if (!_mesh.getProperties().hasPropertyVector("integration_point_offsets"))
            OGS_FATAL(
                "integration_point_data field exists in the input but there is "
                "no integration_point_offsets cell data.");

        auto const& data =
            *_mesh.getProperties().template getPropertyVector<char>(
                "integration_point_data");
        assert(data.getMeshItemType() ==
               MeshLib::MeshItemType::IntegrationPoint);

        auto const& offsets =
            *_mesh.getProperties().template getPropertyVector<std::size_t>(
                "integration_point_offsets");
        assert(offsets.getMeshItemType() == MeshLib::MeshItemType::Cell);

        std::vector<char> local_data;
        assert(_local_assemblers.size() == offsets.size());
        // Starting counting from one; the last cell is handled after the loop.
        std::size_t i = 0;
        for (; i < _local_assemblers.size() - 1; ++i)
        {
            std::size_t const size = offsets[i + 1] - offsets[i];
            local_data.resize(size);
            std::memcpy(local_data.data(), &data[offsets[i]], size);
            _local_assemblers[i]->readIntegrationPointData(local_data);
        }
        {   // last cell
            std::size_t const size = data.size() - offsets[i];
            local_data.resize(size);
            std::memcpy(local_data.data(), &data[offsets[i]], size);
            _local_assemblers[i]->readIntegrationPointData(local_data);
        }
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
