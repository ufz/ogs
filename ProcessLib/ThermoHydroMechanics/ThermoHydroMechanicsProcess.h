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

#include "MeshLib/Elements/Utils.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/HydroMechanics//CreateLocalAssemblers.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "ThermoHydroMechanicsFEM.h"
#include "ThermoHydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoHydroMechanics
{

/**
 * \brief A class to simulate thermo-hydro-mechanical process
 * described by
 *
 * \f[
 *     \mathrm{div} \left[ (\boldsymbol{\mathrm{u}}_\mathrm{S})'_\mathrm{S} + \phi_\mathrm{F} \boldsymbol{\mathrm{w}}_\mathrm{FS} \right]
 *     = \underbrace{\beta^\mathrm{eff}_\mathrm{T} T'_\mathrm{S}}_{\text{first term}} +
 *     \underbrace{\phi_\mathrm{F} \beta_\mathrm{TF} \mathrm{grad}\, T \cdot \boldsymbol{\mathrm{w}}_\mathrm{FS}}_{\text{second term}}
 * \f]
 * \f[
 *    \mathrm{div} \left[ \boldsymbol{\sigma}^\mathrm{E}_\mathrm{S} - \alpha_\mathrm{B} p \boldsymbol{I} \right]
 *    + \varrho^\mathrm{eff} \boldsymbol{g} = \boldsymbol{0}
 * \f]
 * \f[
 *   (\varrho c_p)^\mathrm{eff} \frac{\partial T}{\partial t} + \phi_\mathrm{F} \varrho_\mathrm{FR} c_{p\mathrm{F}} \mathrm{grad}\, T \cdot \boldsymbol{\mathrm{w}}_\mathrm{FS}
 *   - \mathrm{div} \left[ \boldsymbol{\lambda}^\mathrm{eff} \mathrm{grad}\, T \right]
 * \f]
 * where
 *    \f{eqnarray*}{
 *       &\alpha_\mathrm{TS}:&                 \mbox{linear coefficient of thermal expansion of the solid phase,}\\
 *       &\beta_\mathrm{TF}:&                  \mbox{volumetric coefficient of thermal expansion of the fluid phase,}\\
 *       &\varrho:&                            \mbox{density,}\\
 *       &\phi_\mrm{F}:&                       \mbox{porosity,}\\
 *       &\alpha_\mrm{B}:&                     \mbox{Biot coefficient,}\\
 *       &\boldsymbol{\kappa}_\mathrm{F}:&     \mbox{intrinsic permeability}\\
 *       &\mu_\mathrm{FR}:&                    \mbox{viscosity}\\
 *    \f}
 *
 * Detailed model description can refer
 * <a href="Zheng_THMOGS6.pdf" target="_blank"><b>Phase field method</b></a>
 */
template <int DisplacementDim>
class ThermoHydroMechanicsProcess final : public Process
{
    using Base = Process;

public:
    ThermoHydroMechanicsProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        ThermoHydroMechanicsProcessData<DisplacementDim>&& process_data,
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

    void constructDofTable() override
    {
        // Create single component dof in every of the mesh's nodes.
        _mesh_subset_all_nodes.reset(
            new MeshLib::MeshSubset(_mesh, &_mesh.getNodes()));
        // Create single component dof in the mesh's base nodes.
        _base_nodes = MeshLib::getBaseNodes(_mesh.getElements());
        _mesh_subset_base_nodes.reset(
            new MeshLib::MeshSubset(_mesh, &_base_nodes));

        // Collect the mesh subsets in a vector.

        // For temperature, which is the first
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets;
        all_mesh_subsets.emplace_back(_mesh_subset_base_nodes.get());

        // For pressure, which is the second
        all_mesh_subsets.emplace_back(_mesh_subset_base_nodes.get());

        // For displacement.
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            getProcessVariables()[2].get().getNumberOfComponents(),
            [&]() {
                return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
            });

        std::vector<unsigned> const vec_n_components{1, 1, DisplacementDim};
        _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
        std::move(all_mesh_subsets), vec_n_components,
        NumLib::ComponentOrder::BY_LOCATION);
    }

    using LocalAssemblerInterface = ThermoHydroMechanicsLocalAssemblerInterface;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override
    {
        ProcessLib::HydroMechanics::createLocalAssemblers<
            DisplacementDim, ThermoHydroMechanicsLocalAssembler>(
            mesh.getDimension(), mesh.getElements(), dof_table,
            // use displacment process variable for shapefunction order
            getProcessVariables()[2].get().getShapeFunctionOrder(),
            _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
            _process_data);

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
                &ThermoHydroMechanicsLocalAssemblerInterface::getIntPtSigmaXX));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_yy",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoHydroMechanicsLocalAssemblerInterface::getIntPtSigmaYY));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_zz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoHydroMechanicsLocalAssemblerInterface::getIntPtSigmaZZ));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xy",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoHydroMechanicsLocalAssemblerInterface::getIntPtSigmaXY));

        if (DisplacementDim == 3)
        {
            Base::_secondary_variables.addSecondaryVariable(
                "sigma_xz",
                makeExtrapolator(
                    1, getExtrapolator(), _local_assemblers,
                                 &ThermoHydroMechanicsLocalAssemblerInterface::
                                     getIntPtSigmaXZ));

            Base::_secondary_variables.addSecondaryVariable(
                "sigma_yz",
                makeExtrapolator(
                    1, getExtrapolator(), _local_assemblers,
                                 &ThermoHydroMechanicsLocalAssemblerInterface::
                                     getIntPtSigmaYZ));
        }

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xx",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                             &ThermoHydroMechanicsLocalAssemblerInterface::
                                 getIntPtEpsilonXX));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_yy",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                             &ThermoHydroMechanicsLocalAssemblerInterface::
                                 getIntPtEpsilonYY));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_zz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                             &ThermoHydroMechanicsLocalAssemblerInterface::
                                 getIntPtEpsilonZZ));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xy",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                             &ThermoHydroMechanicsLocalAssemblerInterface::
                                 getIntPtEpsilonXY));

        Base::_secondary_variables.addSecondaryVariable(
            "velocity_x",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                             &ThermoHydroMechanicsLocalAssemblerInterface::
                                 getIntPtDarcyVelocityX));

        Base::_secondary_variables.addSecondaryVariable(
            "velocity_y",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                             &ThermoHydroMechanicsLocalAssemblerInterface::
                                 getIntPtDarcyVelocityY));

        Base::_secondary_variables.addSecondaryVariable(
            "velocity_z",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                             &ThermoHydroMechanicsLocalAssemblerInterface::
                                 getIntPtDarcyVelocityZ));
    }

    void assembleConcreteProcess(
        const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
        GlobalVector& b, StaggeredCouplingTerm const& coupling_term) override
    {
        DBUG("Assemble ThermoHydroMechanicsProcess.");

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
        DBUG("AssembleJacobian ThermoHydroMechanicsProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
            _local_assemblers, *_local_to_global_index_map, t, x, xdot,
            dxdot_dx, dx_dx, M, K, b, Jac, coupling_term);
    }

    void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                                    double const dt) override
    {
        DBUG("PreTimestep ThermoHydroMechanicsProcess.");

        _process_data.dt = dt;
        _process_data.t = t;

        GlobalExecutor::executeMemberOnDereferenced(
            &ThermoHydroMechanicsLocalAssemblerInterface::preTimestep,
            _local_assemblers, *_local_to_global_index_map, x, t, dt);
    }

    void postTimestepConcreteProcess(GlobalVector const& x) override
    {
        DBUG("PostTimestep ThermoHydroMechanicsProcess.");

        GlobalExecutor::executeMemberOnDereferenced(
            &ThermoHydroMechanicsLocalAssemblerInterface::postTimestep,
            _local_assemblers, *_local_to_global_index_map, x);
    }

private:
    std::vector<MeshLib::Node*> _base_nodes;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_base_nodes;
    ThermoHydroMechanicsProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<ThermoHydroMechanicsLocalAssemblerInterface>>
        _local_assemblers;
    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;
};

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
