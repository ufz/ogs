/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
             &_assembly_params.fluid_specific_heat_source},
            //! \ogs_file_param_special{prj__processes__process__TES__fluid_specific_isobaric_heat_capacity}
            {"fluid_specific_isobaric_heat_capacity", &_assembly_params.cpG},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_specific_heat_source}
            {"solid_specific_heat_source",
             &_assembly_params.solid_specific_heat_source},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_heat_conductivity}
            {"solid_heat_conductivity", &_assembly_params.solid_heat_cond},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_specific_isobaric_heat_capacity}
            {"solid_specific_isobaric_heat_capacity", &_assembly_params.cpS},
            //! \ogs_file_param_special{prj__processes__process__TES__tortuosity}
            {"tortuosity", &_assembly_params.tortuosity},
            //! \ogs_file_param_special{prj__processes__process__TES__diffusion_coefficient}
            {"diffusion_coefficient",
             &_assembly_params.diffusion_coefficient_component},
            //! \ogs_file_param_special{prj__processes__process__TES__porosity}
            {"porosity", &_assembly_params.poro},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_density_dry}
            {"solid_density_dry", &_assembly_params.rho_SR_dry},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_density_initial}
            {"solid_density_initial", &_assembly_params.initial_solid_density}};

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
            {"characteristic_pressure", &_assembly_params.trafo_p},
            //! \ogs_file_param_special{prj__processes__process__TES__characteristic_temperature}
            {"characteristic_temperature", &_assembly_params.trafo_T},
            //! \ogs_file_param_special{prj__processes__process__TES__characteristic_vapour_mass_fraction}
            {"characteristic_vapour_mass_fraction", &_assembly_params.trafo_x}};

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
        config.getConfigParameterOptional<double>(
            "solid_hydraulic_permeability"))
    {
        DBUG(
            "setting parameter `solid_hydraulic_permeability' to isotropic "
            "value `{:g}'",
            *par);
        const auto dim = mesh.getDimension();
        _assembly_params.solid_perm_tensor =
            Eigen::MatrixXd::Identity(dim, dim) * (*par);
    }

    // reactive system
    _assembly_params.react_sys = Adsorption::AdsorptionReaction::newInstance(
        //! \ogs_file_param{prj__processes__process__TES__reactive_system}
        config.getConfigSubtree("reactive_system"));

    // debug output
    if (auto const param =
            //! \ogs_file_param{prj__processes__process__TES__output_element_matrices}
        config.getConfigParameterOptional<bool>("output_element_matrices"))
    {
        DBUG("output_element_matrices: {:s}", (*param) ? "true" : "false");

        _assembly_params.output_element_matrices = *param;
    }

    // TODO somewhere else
    /*
    if (auto const param =
    //! \ogs_file_param{prj__processes__process__TES__output_global_matrix}
    config.getConfigParameterOptional<bool>("output_global_matrix"))
    {
        DBUG("output_global_matrix: {:s}", (*param) ? "true" : "false");

        this->_process_output.output_global_matrix = *param;
    }
    */
}

void TESProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::createLocalAssemblers<TESLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table, _local_assemblers,
        NumLib::IntegrationOrder{integration_order}, mesh.isAxiallySymmetric(),
        _assembly_params);

    initializeSecondaryVariables();
}

void TESProcess::initializeSecondaryVariables()
{
    // adds a secondary variables to the collection of all secondary variables.
    auto add2nd =
        [&](std::string const& var_name, SecondaryVariableFunctions&& fcts)
    { _secondary_variables.addSecondaryVariable(var_name, std::move(fcts)); };

    // creates an extrapolator
    auto makeEx =
        [&](unsigned const n_components,
            std::vector<double> const& (TESLocalAssemblerInterface::*method)(
                const double /*t*/,
                std::vector<GlobalVector*> const& /*x*/,
                std::vector<
                    NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
                std::vector<double>& /*cache*/)
                const) -> SecondaryVariableFunctions
    {
        return ProcessLib::makeExtrapolator(n_components, getExtrapolator(),
                                            _local_assemblers, method);
    };

    add2nd("solid_density",
           makeEx(1, &TESLocalAssemblerInterface::getIntPtSolidDensity));

    add2nd("reaction_rate",
           makeEx(1, &TESLocalAssemblerInterface::getIntPtReactionRate));

    add2nd("darcy_velocity",
           makeEx(_mesh.getDimension(),
                  &TESLocalAssemblerInterface::getIntPtDarcyVelocity));

    add2nd("loading", makeEx(1, &TESLocalAssemblerInterface::getIntPtLoading));
    add2nd(
        "reaction_damping_factor",
        makeEx(1, &TESLocalAssemblerInterface::getIntPtReactionDampingFactor));

    add2nd("vapour_partial_pressure",
           {1,
            [&](auto&&... args) -> GlobalVector const&
            { return computeVapourPartialPressure(args...); },
            nullptr});
    add2nd("relative_humidity",
           {1,
            [&](auto&&... args) -> GlobalVector const&
            { return computeRelativeHumidity(args...); },
            nullptr});
    add2nd("equilibrium_loading",
           {1,
            [&](auto&&... args) -> GlobalVector const&
            { return computeEquilibriumLoading(args...); },
            nullptr});
}

void TESProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble TESProcess.");

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_table = {
        _local_to_global_index_map.get()};

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        getActiveElementIDs(), dof_table, t, dt, x, x_prev, process_id, &M, &K,
        &b);
}

void TESProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& x_prev, int const process_id,
    GlobalVector& b, GlobalMatrix& Jac)
{
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_table = {
        _local_to_global_index_map.get()};

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, getActiveElementIDs(), dof_table, t, dt, x, x_prev,
        process_id, &b, &Jac);
}

void TESProcess::preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                            const double t,
                                            const double delta_t,
                                            const int process_id)
{
    DBUG("new timestep");

    _assembly_params.delta_t = delta_t;
    _assembly_params.current_time = t;
    ++_assembly_params.timestep;  // TODO remove that

    _x_previous_timestep =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(*x[process_id]);
}

void TESProcess::preIterationConcreteProcess(const unsigned iter,
                                             GlobalVector const& /*x*/)
{
    _assembly_params.iteration_in_current_timestep = iter;
    ++_assembly_params.total_iteration;
    ++_assembly_params.number_of_try_of_iteration;
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

        MathLib::LinAlg::setLocalAccessibleVector(*_x_previous_timestep);

        auto check_variable_bounds =
            [&](std::size_t id, TESLocalAssemblerInterface& loc_asm)
        {
            auto const r_c_indices = NumLib::getRowColumnIndices(
                id, *this->_local_to_global_index_map, indices_cache);
            local_x_cache = x.get(r_c_indices.rows);
            local_x_prev_ts_cache = _x_previous_timestep->get(r_c_indices.rows);

            if (!loc_asm.checkBounds(local_x_cache, local_x_prev_ts_cache))
            {
                check_passed = false;
            }
        };

        GlobalExecutor::executeDereferenced(check_variable_bounds,
                                            _local_assemblers);
    }

    if (!check_passed)
    {
        return NumLib::IterationResult::REPEAT_ITERATION;
    }

    // TODO remove
    DBUG("ts {:d} iteration {:d} (in current ts: {:d}) try {:d} accepted",
         _assembly_params.timestep, _assembly_params.total_iteration,
         _assembly_params.iteration_in_current_timestep,
         _assembly_params.number_of_try_of_iteration);

    _assembly_params.number_of_try_of_iteration = 0;

    return NumLib::IterationResult::SUCCESS;
}

GlobalVector const& TESProcess::computeVapourPartialPressure(
    const double /*t*/,
    std::vector<GlobalVector*> const& x,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
    std::unique_ptr<GlobalVector>& result_cache)
{
    constexpr int process_id = 0;  // monolithic scheme.
    assert(dof_table[process_id] == _local_to_global_index_map.get());

    auto const& dof_table_single = getSingleComponentDOFTable();
    result_cache = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {dof_table_single.dofSizeWithoutGhosts(),
         dof_table_single.dofSizeWithoutGhosts(),
         &dof_table_single.getGhostIndices(), nullptr});

    for (auto const& location :
         MeshLib::views::meshLocations(_mesh, MeshLib::MeshItemType::Node))
    {
        auto const p = NumLib::getNodalValue(*x[process_id], location,
                                             *dof_table[process_id],
                                             COMPONENT_ID_PRESSURE);
        auto const x_mV = NumLib::getNodalValue(*x[process_id], location,
                                                *dof_table[process_id],
                                                COMPONENT_ID_MASS_FRACTION);

        auto const x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, _assembly_params.M_react, _assembly_params.M_inert);

        auto const node_id = location.item_id;
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
    assert(dof_table[process_id] == _local_to_global_index_map.get());

    auto const& dof_table_single = getSingleComponentDOFTable();
    result_cache = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {dof_table_single.dofSizeWithoutGhosts(),
         dof_table_single.dofSizeWithoutGhosts(),
         &dof_table_single.getGhostIndices(), nullptr});

    auto const& x = *xs[0];  // monolithic process
    for (auto const& location :
         MeshLib::views::meshLocations(_mesh, MeshLib::MeshItemType::Node))
    {
        auto const p = NumLib::getNodalValue(
            x, location, *dof_table[process_id], COMPONENT_ID_PRESSURE);
        auto const T = NumLib::getNodalValue(
            x, location, *dof_table[process_id], COMPONENT_ID_TEMPERATURE);
        auto const x_mV = NumLib::getNodalValue(
            x, location, *dof_table[process_id], COMPONENT_ID_MASS_FRACTION);

        auto const x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, _assembly_params.M_react, _assembly_params.M_inert);

        auto const p_S =
            Adsorption::AdsorptionReaction::getEquilibriumVapourPressure(T);

        auto const node_id = location.item_id;
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
    assert(dof_table[process_id] == _local_to_global_index_map.get());

    auto const& dof_table_single = getSingleComponentDOFTable();
    result_cache = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {dof_table_single.dofSizeWithoutGhosts(),
         dof_table_single.dofSizeWithoutGhosts(),
         &dof_table_single.getGhostIndices(), nullptr});

    auto const& x = *xs[0];  // monolithic process
    for (auto const& location :
         MeshLib::views::meshLocations(_mesh, MeshLib::MeshItemType::Node))
    {
        auto const p = NumLib::getNodalValue(
            x, location, *dof_table[process_id], COMPONENT_ID_PRESSURE);
        auto const T = NumLib::getNodalValue(
            x, location, *dof_table[process_id], COMPONENT_ID_TEMPERATURE);
        auto const x_mV = NumLib::getNodalValue(
            x, location, *dof_table[process_id], COMPONENT_ID_MASS_FRACTION);

        auto const x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, _assembly_params.M_react, _assembly_params.M_inert);

        auto const p_V = p * x_nV;
        auto const C_eq =
            (p_V <= 0.0) ? 0.0
                         : _assembly_params.react_sys->getEquilibriumLoading(
                               p_V, T, _assembly_params.M_react);

        auto const node_id = location.item_id;
        result_cache->set(node_id, C_eq);
    }

    return *result_cache;
}

}  // namespace TES

}  // namespace ProcessLib
