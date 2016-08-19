/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TESProcess.h"
#include "BaseLib/Functional.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

// TODO Copied from VectorMatrixAssembler. Could be provided by the DOF table.
inline NumLib::LocalToGlobalIndexMap::RowColumnIndices
getRowColumnIndices_(std::size_t const id,
                     NumLib::LocalToGlobalIndexMap const& dof_table,
                     std::vector<GlobalIndexType>& indices)
{
    assert(dof_table.size() > id);
    indices.clear();

    // Local matrices and vectors will always be ordered by component,
    // no matter what the order of the global matrix is.
    for (unsigned c = 0; c < dof_table.getNumberOfComponents(); ++c)
    {
        auto const& idcs = dof_table(id, c).rows;
        indices.reserve(indices.size() + idcs.size());
        indices.insert(indices.end(), idcs.begin(), idcs.end());
    }

    return NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices,
                                                                 indices);
}

void getVectorValues(
    GlobalVector const& x,
    NumLib::LocalToGlobalIndexMap::RowColumnIndices const& r_c_indices,
    std::vector<double>& local_x)
{
    auto const& indices = r_c_indices.rows;
    local_x.clear();
    local_x.reserve(indices.size());

    for (auto i : indices)
    {
        local_x.emplace_back(x.get(i));
    }
}

namespace ProcessLib
{
namespace TES
{
TESProcess::TESProcess(
    MeshLib::Mesh& mesh,
    Process::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<Process::TimeDiscretization>&& time_discretization,
    std::unique_ptr<NumLib::ConvergenceCriterion>&& convergence_criterion,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    SecondaryVariableCollection&& secondary_variables,
    ProcessOutput&& process_output,
    NumLib::NamedFunctionCaller&& named_function_caller,
    const BaseLib::ConfigTree& config)
    : Process(mesh, nonlinear_solver, std::move(time_discretization),
              std::move(convergence_criterion), std::move(process_variables),
              std::move(secondary_variables), std::move(process_output),
              std::move(named_function_caller))
{
    DBUG("Create TESProcess.");

    // physical parameters for local assembly
    {
        std::vector<std::pair<std::string, double*>> params{
            //! \ogs_file_param_special{process__TES__fluid_specific_heat_source}
            {"fluid_specific_heat_source",
             &_assembly_params.fluid_specific_heat_source},
            //! \ogs_file_param_special{process__TES__fluid_specific_isobaric_heat_capacity}
            {"fluid_specific_isobaric_heat_capacity", &_assembly_params.cpG},
            //! \ogs_file_param_special{process__TES__solid_specific_heat_source}
            {"solid_specific_heat_source",
             &_assembly_params.solid_specific_heat_source},
            //! \ogs_file_param_special{process__TES__solid_heat_conductivity}
            {"solid_heat_conductivity", &_assembly_params.solid_heat_cond},
            //! \ogs_file_param_special{process__TES__solid_specific_isobaric_heat_capacity}
            {"solid_specific_isobaric_heat_capacity", &_assembly_params.cpS},
            //! \ogs_file_param_special{process__TES__tortuosity}
            {"tortuosity", &_assembly_params.tortuosity},
            //! \ogs_file_param_special{process__TES__diffusion_coefficient}
            {"diffusion_coefficient",
             &_assembly_params.diffusion_coefficient_component},
            //! \ogs_file_param_special{process__TES__porosity}
            {"porosity", &_assembly_params.poro},
            //! \ogs_file_param_special{process__TES__solid_density_dry}
            {"solid_density_dry", &_assembly_params.rho_SR_dry},
            //! \ogs_file_param_special{process__TES__solid_density_initial}
            {"solid_density_initial", &_assembly_params.initial_solid_density}};

        for (auto const& p : params)
        {
            if (auto const par =
                    //! \ogs_file_special
                    config.getConfigParameterOptional<double>(p.first))
            {
                DBUG("setting parameter `%s' to value `%g'", p.first.c_str(),
                     *par);
                *p.second = *par;
            }
        }
    }

    // characteristic values of primary variables
    {
        std::vector<std::pair<std::string, Trafo*>> const params{
            //! \ogs_file_param_special{process__TES__characteristic_pressure}
            {"characteristic_pressure", &_assembly_params.trafo_p},
            //! \ogs_file_param_special{process__TES__characteristic_temperature}
            {"characteristic_temperature", &_assembly_params.trafo_T},
            //! \ogs_file_param_special{process__TES__characteristic_vapour_mass_fraction}
            {"characteristic_vapour_mass_fraction", &_assembly_params.trafo_x}};

        for (auto const& p : params)
        {
            if (auto const par =
                    //! \ogs_file_special
                    config.getConfigParameterOptional<double>(p.first))
            {
                INFO("setting parameter `%s' to value `%g'", p.first.c_str(),
                     *par);
                *p.second = Trafo{*par};
            }
        }
    }

    // permeability
    if (auto par =
            //! \ogs_file_param{process__TES__solid_hydraulic_permeability}
            config.getConfigParameterOptional<double>("solid_hydraulic_permeability"))
    {
        DBUG(
            "setting parameter `solid_hydraulic_permeability' to isotropic "
            "value `%g'",
            *par);
        const auto dim = mesh.getDimension();
        _assembly_params.solid_perm_tensor =
            Eigen::MatrixXd::Identity(dim, dim) * (*par);
    }

    // reactive system
    _assembly_params.react_sys = Adsorption::AdsorptionReaction::newInstance(
        //! \ogs_file_param{process__TES__reactive_system}
        config.getConfigSubtree("reactive_system"));

    // debug output
    if (auto const param =
            //! \ogs_file_param{process__TES__output_element_matrices}
            config.getConfigParameterOptional<bool>("output_element_matrices"))
    {
        DBUG("output_element_matrices: %s", (*param) ? "true" : "false");

        _assembly_params.output_element_matrices = *param;
    }

    // TODO somewhere else
    /*
    if (auto const param =
    //! \ogs_file_param{process__TES__output_global_matrix}
    config.getConfigParameterOptional<bool>("output_global_matrix"))
    {
        DBUG("output_global_matrix: %s", (*param) ? "true" : "false");

        this->_process_output.output_global_matrix = *param;
    }
    */
}

void TESProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh, unsigned const integration_order)
{
    ProcessLib::createLocalAssemblers<TESLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table, integration_order,
        _local_assemblers, _assembly_params);

    initializeSecondaryVariables();
}

void TESProcess::initializeSecondaryVariables()
{
    // adds a secondary variables to the collection of all secondary variables.
    auto add2nd = [&](std::string const& var_name, unsigned const n_components,
                      SecondaryVariableFunctions&& fcts) {
        _secondary_variables.addSecondaryVariable(var_name, n_components,
                                                  std::move(fcts));
    };

    // creates an extrapolator
    auto makeEx =
        [&](std::vector<double> const& (TESLocalAssemblerInterface::*method)(
            std::vector<double>&)const) -> SecondaryVariableFunctions {
        return ProcessLib::makeExtrapolator(getExtrapolator(),
                                            _local_assemblers, method);
    };

    // named functions: vapour partial pressure ////////////////////////////////
    auto p_V_fct = [=](const double p, const double x_mV) {
        const double x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, _assembly_params.M_react, _assembly_params.M_inert);
        return p*x_nV;
    };
    _named_function_caller.addNamedFunction(
        {"vapour_partial_pressure",
         {"pressure", "vapour_mass_fraction"},
         BaseLib::easyBind(std::move(p_V_fct))});
    _named_function_caller.plug("vapour_partial_pressure", "pressure",
                                "TES_pressure");
    _named_function_caller.plug("vapour_partial_pressure",
                                "vapour_mass_fraction",
                                "TES_vapour_mass_fraction");
    // /////////////////////////////////////////////////////////////////////////

    // named functions: solid density //////////////////////////////////////////
    std::unique_ptr<CachedSecondaryVariable> solid_density(
        new CachedSecondaryVariable(
            "TES_solid_density", getExtrapolator(), _local_assemblers,
            &TESLocalAssemblerInterface::getIntPtSolidDensity,
            _secondary_variable_context));

    for (auto&& fct : solid_density->getNamedFunctions())
        _named_function_caller.addNamedFunction(std::move(fct));

    add2nd("solid_density", 1, solid_density->getExtrapolator());

    _cached_secondary_variables.emplace_back(std::move(solid_density));
    // /////////////////////////////////////////////////////////////////////////

    add2nd("reaction_rate", 1,
           makeEx(&TESLocalAssemblerInterface::getIntPtReactionRate));

    add2nd("velocity_x", 1,
           makeEx(&TESLocalAssemblerInterface::getIntPtDarcyVelocityX));
    if (_mesh.getDimension() >= 2)
        add2nd("velocity_y", 1,
               makeEx(&TESLocalAssemblerInterface::getIntPtDarcyVelocityY));
    if (_mesh.getDimension() >= 3)
        add2nd("velocity_z", 1,
               makeEx(&TESLocalAssemblerInterface::getIntPtDarcyVelocityZ));

    add2nd("loading", 1, makeEx(&TESLocalAssemblerInterface::getIntPtLoading));
    add2nd("reaction_damping_factor", 1,
          makeEx(&TESLocalAssemblerInterface::getIntPtReactionDampingFactor));

    add2nd("relative_humidity", 1,
           {BaseLib::easyBind(&TESProcess::computeRelativeHumidity, this),
            nullptr});
    add2nd("equilibrium_loading", 1,
           {BaseLib::easyBind(&TESProcess::computeEquilibriumLoading, this),
            nullptr});
}

void TESProcess::assembleConcreteProcess(const double t,
                                                      GlobalVector const& x,
                                                      GlobalMatrix& M,
                                                      GlobalMatrix& K,
                                                      GlobalVector& b)
{
    DBUG("Assemble TESProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberOnDereferenced(
        &TESLocalAssemblerInterface::assemble, _local_assemblers,
        *_local_to_global_index_map, t, x, M, K, b);
}

void TESProcess::preTimestep(GlobalVector const& x, const double t,
                                          const double delta_t)
{
    DBUG("new timestep");

    _assembly_params.delta_t = delta_t;
    _assembly_params.current_time = t;
    ++_assembly_params.timestep;  // TODO remove that

    _x_previous_timestep =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(x);
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
    if (this->_process_output.output_iteration_results)
    {
        DBUG("output results of iteration %li",
             _assembly_params.total_iteration);
        std::string fn =
            "tes_iter_" + std::to_string(_assembly_params.total_iteration) +
            +"_ts_" + std::to_string(_assembly_params.timestep) + "_" +
            std::to_string(_assembly_params.iteration_in_current_timestep) +
            "_" + std::to_string(_assembly_params.number_of_try_of_iteration) +
            ".vtu";

        this->output(fn, 0, x);
    }

    bool check_passed = true;

    if (!Trafo::constrained)
    {
        // bounds checking only has to happen if the vapour mass fraction is
        // non-logarithmic.

        std::vector<GlobalIndexType> indices_cache;
        std::vector<double> local_x_cache;
        std::vector<double> local_x_prev_ts_cache;

        auto check_variable_bounds = [&](std::size_t id,
                                         TESLocalAssemblerInterface& loc_asm) {
            auto const r_c_indices = getRowColumnIndices_(
                id, *this->_local_to_global_index_map, indices_cache);
            getVectorValues(x, r_c_indices, local_x_cache);
            getVectorValues(*_x_previous_timestep, r_c_indices,
                            local_x_prev_ts_cache);

            if (!loc_asm.checkBounds(local_x_cache, local_x_prev_ts_cache))
                check_passed = false;
        };

        GlobalExecutor::executeDereferenced(check_variable_bounds,
                                         _local_assemblers);
    }

    if (!check_passed)
        return NumLib::IterationResult::REPEAT_ITERATION;

    // TODO remove
    DBUG("ts %lu iteration %lu (in current ts: %lu) try %u accepted",
         _assembly_params.timestep, _assembly_params.total_iteration,
         _assembly_params.iteration_in_current_timestep,
         _assembly_params.number_of_try_of_iteration);

    _assembly_params.number_of_try_of_iteration = 0;

    return NumLib::IterationResult::SUCCESS;
}

GlobalVector const&
TESProcess::computeVapourPartialPressure(
    GlobalVector const& x,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::unique_ptr<GlobalVector>& result_cache)
{
    assert(&dof_table == _local_to_global_index_map.get());

    auto const& dof_table_single = getSingleComponentDOFTable();
    result_cache = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {dof_table_single.dofSizeWithoutGhosts(),
         dof_table_single.dofSizeWithoutGhosts(),
         &dof_table_single.getGhostIndices(), nullptr});

    GlobalIndexType const nnodes = _mesh.getNumberOfNodes();

    for (GlobalIndexType node_id = 0; node_id < nnodes; ++node_id)
    {
        auto const p = NumLib::getNodalValue(x, _mesh, dof_table, node_id,
                                             COMPONENT_ID_PRESSURE);
        auto const x_mV = NumLib::getNodalValue(x, _mesh, dof_table, node_id,
                                                COMPONENT_ID_MASS_FRACTION);

        auto const x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, _assembly_params.M_react, _assembly_params.M_inert);

        // TODO Problems with PETSc? (local vs. global index)
        result_cache->set(node_id, p * x_nV);
    }

    return *result_cache;
}

GlobalVector const&
TESProcess::computeRelativeHumidity(
    GlobalVector const& x,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::unique_ptr<GlobalVector>& result_cache)
{
    assert(&dof_table == _local_to_global_index_map.get());

    auto const& dof_table_single = getSingleComponentDOFTable();
    result_cache = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {dof_table_single.dofSizeWithoutGhosts(),
         dof_table_single.dofSizeWithoutGhosts(),
         &dof_table_single.getGhostIndices(), nullptr});

    GlobalIndexType const nnodes = _mesh.getNumberOfNodes();

    for (GlobalIndexType node_id = 0; node_id < nnodes; ++node_id)
    {
        auto const p = NumLib::getNodalValue(x, _mesh, dof_table, node_id,
                                             COMPONENT_ID_PRESSURE);
        auto const T = NumLib::getNodalValue(x, _mesh, dof_table, node_id,
                                             COMPONENT_ID_TEMPERATURE);
        auto const x_mV = NumLib::getNodalValue(x, _mesh, dof_table, node_id,
                                                COMPONENT_ID_MASS_FRACTION);

        auto const x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, _assembly_params.M_react, _assembly_params.M_inert);

        auto const p_S =
            Adsorption::AdsorptionReaction::getEquilibriumVapourPressure(T);

        // TODO Problems with PETSc? (local vs. global index)
        result_cache->set(node_id, p * x_nV / p_S);
    }

    return *result_cache;
}

GlobalVector const&
TESProcess::computeEquilibriumLoading(
    GlobalVector const& x,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::unique_ptr<GlobalVector>& result_cache)
{
    assert(&dof_table == _local_to_global_index_map.get());

    auto const& dof_table_single = getSingleComponentDOFTable();
    result_cache = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {dof_table_single.dofSizeWithoutGhosts(),
         dof_table_single.dofSizeWithoutGhosts(),
         &dof_table_single.getGhostIndices(), nullptr});

    GlobalIndexType const nnodes = _mesh.getNumberOfNodes();

    for (GlobalIndexType node_id = 0; node_id < nnodes; ++node_id)
    {
        auto const p = NumLib::getNodalValue(x, _mesh, dof_table, node_id,
                                             COMPONENT_ID_PRESSURE);
        auto const T = NumLib::getNodalValue(x, _mesh, dof_table, node_id,
                                             COMPONENT_ID_TEMPERATURE);
        auto const x_mV = NumLib::getNodalValue(
            x, _mesh, dof_table, node_id, COMPONENT_ID_MASS_FRACTION);

        auto const x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, _assembly_params.M_react, _assembly_params.M_inert);

        auto const p_V = p * x_nV;
        auto const C_eq =
            (p_V <= 0.0) ? 0.0
                         : _assembly_params.react_sys->getEquilibriumLoading(
                               p_V, T, _assembly_params.M_react);

        // TODO Problems with PETSc? (local vs. global index)
        result_cache->set(node_id, C_eq);
    }

    return *result_cache;
}

}  // namespace TES

}  // namespace ProcessLib
