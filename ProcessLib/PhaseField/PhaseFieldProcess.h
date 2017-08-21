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

#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

#include "PhaseFieldFEM.h"
#include "PhaseFieldProcessData.h"

namespace ProcessLib
{
namespace PhaseField
{

/**
 * \brief A class to simulate mechanical fracturing process
 * using phase-field approach in solids described by
 *
 * \f[
 *     \mathrm{div} \left[ \left(d^2 + k \right) \boldsymbol{\sigma}_0^+
 *      + \boldsymbol{\sigma}_0^- \right] +
 *        \varrho \boldsymbol{b} = \boldsymbol{0}
 * \f]
 * \f[
 *     2d H(\boldsymbol{\epsilon}_\mathrm{el})
 *      - \frac{1 - d}{2 \varepsilon} g_\mathrm{c}
 *      - 2 \varepsilon g_\mathrm{c} \mathrm{div}(\mathrm{grad} d) = 0
 * \f]
 * where
 *    \f{eqnarray*}{
 *       &d:&                    \mbox{order parameter,}\\
 *       &\varrho:&              \mbox{density,}\\
 *       &H:&                    \mbox{history field,}\\
 *       &g_\mathrm{c}:&         \mbox{fracture energy,}\\
 *       &\varepsilon:&          \mbox{length scale}\\
 *    \f}
 *
 * Detailed model description can refer
 * <a href="Miao_Biot2017.pdff" target="_blank"><b>Phase field method</b></a>
 */
template <int DisplacementDim>
class PhaseFieldProcess final : public Process
{
    using Base = Process;

public:
    PhaseFieldProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        PhaseFieldProcessData<DisplacementDim>&& process_data,
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
    using LocalAssemblerInterface = PhaseFieldLocalAssemblerInterface;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override
    {
        ProcessLib::SmallDeformation::createLocalAssemblers<
                DisplacementDim, PhaseFieldLocalAssembler>(
                mesh.getElements(), dof_table, _local_assemblers,
                mesh.isAxiallySymmetric(), integration_order, _process_data);

        // TODO move the two data members somewhere else.
        // for extrapolation of secondary variables
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets_single_component;
        all_mesh_subsets_single_component.emplace_back(
            _mesh_subset_all_nodes.get());
        _local_to_global_index_map_single_component.reset(
            new NumLib::LocalToGlobalIndexMap(
                std::move(all_mesh_subsets_single_component),
                // by location order is needed for output
                NumLib::ComponentOrder::BY_LOCATION));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xx",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &PhaseFieldLocalAssemblerInterface::getIntPtSigmaXX));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_yy",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &PhaseFieldLocalAssemblerInterface::getIntPtSigmaYY));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_zz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &PhaseFieldLocalAssemblerInterface::getIntPtSigmaZZ));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xy",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &PhaseFieldLocalAssemblerInterface::getIntPtSigmaXY));

        if (DisplacementDim == 3)
        {
            Base::_secondary_variables.addSecondaryVariable(
                "sigma_xz",
                makeExtrapolator(
                    1, getExtrapolator(), _local_assemblers,
                    &PhaseFieldLocalAssemblerInterface::getIntPtSigmaXZ));

            Base::_secondary_variables.addSecondaryVariable(
                "sigma_yz",
                makeExtrapolator(
                    1, getExtrapolator(), _local_assemblers,
                    &PhaseFieldLocalAssemblerInterface::getIntPtSigmaYZ));
        }
        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xx",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &PhaseFieldLocalAssemblerInterface::getIntPtEpsilonXX));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_yy",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &PhaseFieldLocalAssemblerInterface::getIntPtEpsilonYY));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_zz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &PhaseFieldLocalAssemblerInterface::getIntPtEpsilonZZ));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xy",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &PhaseFieldLocalAssemblerInterface::getIntPtEpsilonXY));
        if (DisplacementDim == 3)
        {
            Base::_secondary_variables.addSecondaryVariable(
                "epsilon_yz",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                                 &PhaseFieldLocalAssemblerInterface::
                                     getIntPtEpsilonYZ));

            Base::_secondary_variables.addSecondaryVariable(
                "epsilon_xz",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                                 &PhaseFieldLocalAssemblerInterface::
                                     getIntPtEpsilonXZ));
        }
    }

    void assembleConcreteProcess(
        const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
        GlobalVector& b) override
    {
        DBUG("Assemble PhaseFieldProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assemble,
            _local_assemblers, *_local_to_global_index_map, t, x, M, K, b,
            *_coupling_term.lock());
    }

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override
    {
        // DBUG("AssembleJacobian PhaseFieldProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
            _local_assemblers, *_local_to_global_index_map, t, x, xdot,
            dxdot_dx, dx_dx, M, K, b, Jac, *_coupling_term.lock());
    }

    void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                                    double const dt) override
    {
        DBUG("PreTimestep PhaseFieldProcess.");

        _process_data.dt = dt;
        _process_data.t = t;

        GlobalExecutor::executeMemberOnDereferenced(
            &PhaseFieldLocalAssemblerInterface::preTimestep, _local_assemblers,
            *_local_to_global_index_map, x, t, dt);
    }

    void postTimestepConcreteProcess(GlobalVector const& x) override
    {
        DBUG("PostTimestep PhaseFieldProcess.");

        GlobalExecutor::executeMemberOnDereferenced(
            &PhaseFieldLocalAssemblerInterface::postTimestep, _local_assemblers,
            *_local_to_global_index_map, x);
    }

private:
    PhaseFieldProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;
};

}  // namespace PhaseField
}  // namespace ProcessLib
