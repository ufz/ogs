/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TESProcess.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace TES
{
TESProcess::TESProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    SecondaryVariableCollection&& secondary_variables,
    const BaseLib::ConfigTree& config)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables))
{
    DBUG("Create TESProcess.");

    // physical parameters for local assembly
    {
        std::vector<std::pair<std::string, double*>> params{
            //! \ogs_file_param_special{prj__processes__process__TES__fluid_specific_heat_source}
            {"fluid_specific_heat_source",
             &assembly_params_.fluid_specific_heat_source},
            //! \ogs_file_param_special{prj__processes__process__TES__fluid_specific_isobaric_heat_capacity}
            {"fluid_specific_isobaric_heat_capacity", &assembly_params_.cpG},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_specific_heat_source}
            {"solid_specific_heat_source",
             &assembly_params_.solid_specific_heat_source},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_heat_conductivity}
            {"solid_heat_conductivity", &assembly_params_.solid_heat_cond},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_specific_isobaric_heat_capacity}
            {"solid_specific_isobaric_heat_capacity", &assembly_params_.cpS},
            //! \ogs_file_param_special{prj__processes__process__TES__tortuosity}
            {"tortuosity", &assembly_params_.tortuosity},
            //! \ogs_file_param_special{prj__processes__process__TES__diffusion_coefficient}
            {"diffusion_coefficient",
             &assembly_params_.diffusion_coefficient_component},
            //! \ogs_file_param_special{prj__processes__process__TES__porosity}
            {"porosity", &assembly_params_.poro},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_density_dry}
            {"solid_density_dry", &assembly_params_.rho_SR_dry},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_density_initial}
            {"solid_density_initial", &assembly_params_.initial_solid_density}};

        for (auto const& p : params)
        {
            if (auto const par =
                    //! \ogs_file_special
                    config.getConfigParameterOptional<double>(p.first))
            {
                DBUG("setting parameter `{:s}' to value `{:g}'", p.first, *par);
                *p.second = *par;
            }
        }
    }

    // characteristic values of primary variables
    {
        std::vector<std::pair<std::string, Trafo*>> const params{
            //! \ogs_file_param_special{prj__processes__process__TES__characteristic_pressure}
            {"characteristic_pressure", &assembly_params_.trafo_p},
            //! \ogs_file_param_special{prj__processes__process__TES__characteristic_temperature}
            {"characteristic_temperature", &assembly_params_.trafo_T},
            //! \ogs_file_param_special{prj__processes__process__TES__characteristic_vapour_mass_fraction}
            {"characteristic_vapour_mass_fraction", &assembly_params_.trafo_x}};

        for (auto const& p : params)
        {
            if (auto const par =
                    //! \ogs_file_special
                    config.getConfigParameterOptional<double>(p.first))
            {
                INFO("setting parameter `{:s}' to value `{:g}'", p.first, *par);
                *p.second = Trafo{*par};
            }
        }
    }

    // permeability
    if (auto par =
            //! \ogs_file_param{prj__processes__process__TES__solid_hydraulic_permeability}
            config.getConfigParameterOptional<double>("solid_hydraulic_permeability"))
    {
        DBUG(
            "setting parameter `solid_hydraulic_permeability' to isotropic "
            "value `{:g}'",
            *par);
        const auto dim = mesh.getDimension();
        assembly_params_.solid_perm_tensor =
            Eigen::MatrixXd::Identity(dim, dim) * (*par);
    }

    // reactive system
    assembly_params_.react_sys = Adsorption::AdsorptionReaction::newInstance(
        //! \ogs_file_param{prj__processes__process__TES__reactive_system}
        config.getConfigSubtree("reactive_system"));

    // debug output
    if (auto const param =
            //! \ogs_file_param{prj__processes__process__TES__output_element_matrices}
            config.getConfigParameterOptional<bool>("output_element_matrices"))
    {
        DBUG("output_element_matrices: {:s}", (*param) ? "true" : "false");

        assembly_params_.output_element_matrices = *param;
    }

    // TODO somewhere else
    /*
    if (auto const param =
    //! \ogs_file_param{prj__processes__process__TES__output_global_matrix}
    config.getConfigParameterOptional<bool>("output_global_matrix"))
    {
        DBUG("output_global_matrix: {:s}", (*param) ? "true" : "false");

        this->process_output_.output_global_matrix = *param;
    }
    */
}

void TESProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    const int monolithic_process_id = 0;
    ProcessLib::ProcessVariable const& pv =
        getProcessVariables(monolithic_process_id)[0];
    ProcessLib::createLocalAssemblers<TESLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), local_assemblers_,
        mesh.isAxiallySymmetric(), integration_order, assembly_params_);

    initializeSecondaryVariables();
}

void TESProcess::initializeSecondaryVariables()
{
    // adds a secondary variables to the collection of all secondary variables.
    auto add2nd = [&](std::string const& var_name,
                      SecondaryVariableFunctions&& fcts) {
        secondary_variables_.addSecondaryVariable(var_name, std::move(fcts));
    };

    // creates an extrapolator
    auto makeEx =
        [&](unsigned const n_components,
            std::vector<double> const& (TESLocalAssemblerInterface::*method)(
                const double /*t*/,
                std::vector<GlobalVector*> const& /*x*/,
                std::vector<
                    NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
                std::vector<double>& /*cache*/)
                const) -> SecondaryVariableFunctions {
        return ProcessLib::makeExtrapolator(n_components, getExtrapolator(),
                                            local_assemblers_, method);
    };

    add2nd("solid_density",
           makeEx(1, &TESLocalAssemblerInterface::getIntPtSolidDensity));

    add2nd("reaction_rate",
           makeEx(1, &TESLocalAssemblerInterface::getIntPtReactionRate));

    add2nd("darcy_velocity",
           makeEx(mesh_.getDimension(),
                  &TESLocalAssemblerInterface::getIntPtDarcyVelocity));

    add2nd("loading", makeEx(1, &TESLocalAssemblerInterface::getIntPtLoading));
    add2nd(
        "reaction_damping_factor",
        makeEx(1, &TESLocalAssemblerInterface::getIntPtReactionDampingFactor));

    add2nd("vapour_partial_pressure",
           {1,
            [&](auto&&... args) -> GlobalVector const& {
                return computeVapourPartialPressure(args...);
            },
            nullptr});
    add2nd("relative_humidity",
           {1,
            [&](auto&&... args) -> GlobalVector const& {
                return computeRelativeHumidity(args...);
            },
            nullptr});
    add2nd("equilibrium_loading",
           {1,
            [&](auto&&... args) -> GlobalVector const& {
                return computeEquilibriumLoading(args...);
            },
            nullptr});
}

void TESProcess::assembleConcreteProcess(const double t, double const dt,
                                         std::vector<GlobalVector*> const& x,
                                         std::vector<GlobalVector*> const& xdot,
                                         int const process_id, GlobalMatrix& M,
                                         GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble TESProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
       dof_table = {std::ref(*local_to_global_index_map_)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assemble, local_assemblers_,
        pv.getActiveElementIDs(), dof_table, t, dt, x, xdot, process_id, M, K,
        b, coupled_solutions_);
}

void TESProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
       dof_table = {std::ref(*local_to_global_index_map_)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        global_assembler_, &VectorMatrixAssembler::assembleWithJacobian,
        local_assemblers_, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, coupled_solutions_);
}

void TESProcess::preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                            const double t,
                                            const double delta_t,
                                            const int process_id)
{
    DBUG("new timestep");

    assembly_params_.delta_t = delta_t;
    assembly_params_.current_time = t;
    ++assembly_params_.timestep;  // TODO remove that

    x_previous_timestep_ =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(*x[process_id]);
}

void TESProcess::preIterationConcreteProcess(const unsigned iter,
                                             GlobalVector const& /*x*/)
{
    assembly_params_.iteration_in_current_timestep = iter;
    ++assembly_params_.total_iteration;
    ++assembly_params_.number_of_try_of_iteration;
}

NumLib::IterationResult TESProcess::postIterationConcreteProcess(
    GlobalVector const& x)
{
    bool check_passed = true;

    if (!Trafo::constrained)
    {
        // bounds checking only has to happen if the vapour mass fraction is
        // non-logarithmic.

        std::vector<GlobalIndexType> indices_cache;
        std::vector<double> local_x_cache;
        std::vector<double> local_x_prev_ts_cache;

        MathLib::LinAlg::setLocalAccessibleVector(*x_previous_timestep_);

        auto check_variable_bounds = [&](std::size_t id,
                                         TESLocalAssemblerInterface& loc_asm) {
            auto const r_c_indices = NumLib::getRowColumnIndices(
                id, *this->local_to_global_index_map_, indices_cache);
            local_x_cache = x.get(r_c_indices.rows);
            local_x_prev_ts_cache = x_previous_timestep_->get(r_c_indices.rows);

            if (!loc_asm.checkBounds(local_x_cache, local_x_prev_ts_cache))
            {
                check_passed = false;
            }
        };

        GlobalExecutor::executeDereferenced(check_variable_bounds,
                                         local_assemblers_);
    }

    if (!check_passed)
    {
        return NumLib::IterationResult::REPEAT_ITERATION;
    }

    // TODO remove
    DBUG("ts {:d} iteration {:d} (in current ts: {:d}) try {:d} accepted",
         assembly_params_.timestep, assembly_params_.total_iteration,
         assembly_params_.iteration_in_current_timestep,
         assembly_params_.number_of_try_of_iteration);

    assembly_params_.number_of_try_of_iteration = 0;

    return NumLib::IterationResult::SUCCESS;
}

GlobalVector const& TESProcess::computeVapourPartialPressure(
    const double /*t*/,
    std::vector<GlobalVector*> const& x,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
    std::unique_ptr<GlobalVector>& result_cache)
{
    constexpr int process_id = 0;  // monolithic scheme.
    assert(dof_table[process_id] == local_to_global_index_map_.get());

    auto const& dof_table_single = getSingleComponentDOFTable();
    result_cache = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {dof_table_single.dofSizeWithoutGhosts(),
         dof_table_single.dofSizeWithoutGhosts(),
         &dof_table_single.getGhostIndices(), nullptr});

    GlobalIndexType const nnodes = mesh_.getNumberOfNodes();

    for (GlobalIndexType node_id = 0; node_id < nnodes; ++node_id)
    {
        auto const p =
            NumLib::getNodalValue(*x[process_id], mesh_, *dof_table[process_id],
                                  node_id, COMPONENT_ID_PRESSURE);
        auto const x_mV =
            NumLib::getNodalValue(*x[process_id], mesh_, *dof_table[process_id],
                                  node_id, COMPONENT_ID_MASS_FRACTION);

        auto const x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, assembly_params_.M_react, assembly_params_.M_inert);

        result_cache->set(node_id, p * x_nV);
    }

    return *result_cache;
}

GlobalVector const& TESProcess::computeRelativeHumidity(
    double const /*t*/,
    std::vector<GlobalVector*> const& xs,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
    std::unique_ptr<GlobalVector>& result_cache)
{
    constexpr int process_id = 0;  // monolithic scheme.
    assert(dof_table[process_id] == local_to_global_index_map_.get());

    auto const& dof_table_single = getSingleComponentDOFTable();
    result_cache = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {dof_table_single.dofSizeWithoutGhosts(),
         dof_table_single.dofSizeWithoutGhosts(),
         &dof_table_single.getGhostIndices(), nullptr});

    GlobalIndexType const nnodes = mesh_.getNumberOfNodes();

    auto const& x = *xs[0];  // monolithic process
    for (GlobalIndexType node_id = 0; node_id < nnodes; ++node_id)
    {
        auto const p = NumLib::getNodalValue(x, mesh_, *dof_table[process_id],
                                             node_id, COMPONENT_ID_PRESSURE);
        auto const T = NumLib::getNodalValue(x, mesh_, *dof_table[process_id],
                                             node_id, COMPONENT_ID_TEMPERATURE);
        auto const x_mV =
            NumLib::getNodalValue(x, mesh_, *dof_table[process_id], node_id,
                                  COMPONENT_ID_MASS_FRACTION);

        auto const x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, assembly_params_.M_react, assembly_params_.M_inert);

        auto const p_S =
            Adsorption::AdsorptionReaction::getEquilibriumVapourPressure(T);

        result_cache->set(node_id, p * x_nV / p_S);
    }

    return *result_cache;
}

GlobalVector const& TESProcess::computeEquilibriumLoading(
    double const /*t*/,
    std::vector<GlobalVector*> const& xs,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
    std::unique_ptr<GlobalVector>& result_cache)
{
    constexpr int process_id = 0;  // monolithic scheme.
    assert(dof_table[process_id] == local_to_global_index_map_.get());

    auto const& dof_table_single = getSingleComponentDOFTable();
    result_cache = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {dof_table_single.dofSizeWithoutGhosts(),
         dof_table_single.dofSizeWithoutGhosts(),
         &dof_table_single.getGhostIndices(), nullptr});

    GlobalIndexType const nnodes = mesh_.getNumberOfNodes();

    auto const& x = *xs[0];  // monolithic process
    for (GlobalIndexType node_id = 0; node_id < nnodes; ++node_id)
    {
        auto const p = NumLib::getNodalValue(x, mesh_, *dof_table[process_id],
                                             node_id, COMPONENT_ID_PRESSURE);
        auto const T = NumLib::getNodalValue(x, mesh_, *dof_table[process_id],
                                             node_id, COMPONENT_ID_TEMPERATURE);
        auto const x_mV =
            NumLib::getNodalValue(x, mesh_, *dof_table[process_id], node_id,
                                  COMPONENT_ID_MASS_FRACTION);

        auto const x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, assembly_params_.M_react, assembly_params_.M_inert);

        auto const p_V = p * x_nV;
        auto const C_eq =
            (p_V <= 0.0) ? 0.0
                         : assembly_params_.react_sys->getEquilibriumLoading(
                               p_V, T, assembly_params_.M_react);

        result_cache->set(node_id, C_eq);
    }

    return *result_cache;
}

}  // namespace TES

}  // namespace ProcessLib
