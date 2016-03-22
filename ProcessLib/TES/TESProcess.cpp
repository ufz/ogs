/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TESProcess.h"

namespace ProcessLib
{

namespace TES
{

template<typename GlobalSetup>
TESProcess<GlobalSetup>::
TESProcess(MeshLib::Mesh& mesh,
           typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
           std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&& time_discretization,
           std::vector<ProcessVariable> const& variables,
           std::vector<std::unique_ptr<ParameterBase>> const& /*parameters*/,
           const BaseLib::ConfigTree& config)
    : Process<GlobalSetup>(mesh, nonlinear_solver, std::move(time_discretization))
{
    DBUG("Create TESProcess.");

    // primary variables
    {
        auto const proc_vars = config.getConfSubtree("process_variables");

        for (auto const& var_name
             : { "fluid_pressure", "temperature", "vapour_mass_fraction" })
        {
            auto& variable = findProcessVariable(proc_vars, var_name, variables);

            // TODO extend Process constructor to cover that
            BP::_process_variables.emplace_back(variable);
        }
    }

    // secondary variables
    if (auto proc_vars = config.getConfSubtreeOptional("secondary_variables"))
    {
        auto add_secondary_variable =
                [this, &proc_vars](
                std::string const& var, const unsigned num_components,
                typename ProcessLib::SecondaryVariableFunctions<GlobalVector>&& fcts)
                -> void
        {
            if (auto variable = proc_vars->getConfParamOptional<std::string>(var))
            {
                BP::_process_output.secondary_variables.push_back(
                    {var, num_components, std::move(fcts)});
            }
        };

        add_secondary_variable("solid_density",  1, makeExtrapolator(TESIntPtVariables::SOLID_DENSITY));
        add_secondary_variable("reaction_rate",  1, makeExtrapolator(TESIntPtVariables::REACTION_RATE));
        add_secondary_variable("velocity_x",     1, makeExtrapolator(TESIntPtVariables::VELOCITY_X));
        if (BP::_mesh.getDimension() >= 2)
            add_secondary_variable("velocity_y", 1, makeExtrapolator(TESIntPtVariables::VELOCITY_Y));
        if (BP::_mesh.getDimension() >= 3)
            add_secondary_variable("velocity_z", 1, makeExtrapolator(TESIntPtVariables::VELOCITY_Z));

        add_secondary_variable("loading",        1, makeExtrapolator(TESIntPtVariables::LOADING));
        add_secondary_variable("reaction_damping_factor", 1,
                               makeExtrapolator(TESIntPtVariables::REACTION_DAMPING_FACTOR));

        namespace PH = std::placeholders;
        using Self = TESProcess<GlobalSetup>;

        add_secondary_variable("vapour_partial_pressure", 1,
            {std::bind(&Self::computeVapourPartialPressure, this, PH::_1, PH::_2), nullptr});
        add_secondary_variable("relative_humidity",       1,
            {std::bind(&Self::computeRelativeHumidity,      this, PH::_1, PH::_2), nullptr});
        add_secondary_variable("equilibrium_loading",     1,
            {std::bind(&Self::computeEquilibriumLoading,    this, PH::_1, PH::_2), nullptr});
    }

    // variables for output
    if (auto output = config.getConfSubtreeOptional("output")) {
        auto& output_variables = BP::_process_output.output_variables;
        if (auto out_vars = output->getConfSubtreeOptional("variables"))
        {
            for (auto out_var : out_vars->getConfParamList<std::string>("variable"))
            {
                if (output_variables.find(out_var) != output_variables.cend())
                {
                    ERR("output variable `%s' specified twice.", out_var.c_str());
                    std::abort();
                }

                auto pred = [&out_var](ProcessVariable const& pv) {
                    return pv.getName() == out_var;
                };

                // check if process variable
                auto const& pcs_var = std::find_if(
                    BP::_process_variables.cbegin(), BP::_process_variables.cend(),
                    pred);

                if (pcs_var == BP::_process_variables.cend())
                {
                    auto pred2 = [&out_var](SecondaryVariable<GlobalVector> const& p) {
                        return p.name == out_var;
                    };

                    // check if secondary variable
                    auto const& pcs_var2 = std::find_if(
                        BP::_process_output.secondary_variables.cbegin(),
                        BP::_process_output.secondary_variables.cend(),
                        pred2);

                    if (pcs_var2 == BP::_process_output.secondary_variables.cend())
                    {
                        ERR("Output variable `%s' is neither a process variable nor a"
                            " secondary variable", out_var.c_str());
                        std::abort();
                    }
                }

                DBUG("adding output variable `%s'", out_var.c_str());
                output_variables.insert(out_var);
            }

            if (auto out_resid = output->getConfParamOptional<bool>("output_extrapolation_residuals"))
            {
                BP::_process_output.output_residuals = *out_resid;
            }
        }
    }

    {
        std::vector<std::pair<const std::string, double*> > params{
            { "fluid_specific_heat_source",            &_assembly_params.fluid_specific_heat_source },
            { "fluid_specific_isobaric_heat_capacity", &_assembly_params.cpG },
            { "solid_specific_heat_source",            &_assembly_params.solid_specific_heat_source },
            { "solid_heat_conductivity",               &_assembly_params.solid_heat_cond },
            { "solid_specific_isobaric_heat_capacity", &_assembly_params.cpS },
            { "tortuosity",                            &_assembly_params.tortuosity },
            { "diffusion_coefficient",                 &_assembly_params.diffusion_coefficient_component },
            { "porosity",                              &_assembly_params.poro },
            { "solid_density_dry",                     &_assembly_params.rho_SR_dry },
            { "solid_density_initial",                 &_assembly_params.initial_solid_density }
        };

        for (auto const& p : params)
        {
            if (auto const par = config.getConfParamOptional<double>(p.first)) {
                DBUG("setting parameter `%s' to value `%g'", p.first.c_str(), *par);
                *p.second = *par;
            }
        }
    }

    // characteristic values of primary variables
    {
        std::vector<std::pair<const std::string, Trafo*> > const params{
            { "characteristic_pressure",             &_assembly_params.trafo_p },
            { "characteristic_temperature",          &_assembly_params.trafo_T },
            { "characteristic_vapour_mass_fraction", &_assembly_params.trafo_x }
        };

        for (auto const& p : params)
        {
            if (auto const par = config.getConfParamOptional<double>(p.first)) {
                INFO("setting parameter `%s' to value `%g'", p.first.c_str(), *par);
                *p.second = Trafo{*par};
            }
        }
    }

    // permeability
    if (auto par = config.getConfParamOptional<double>("solid_hydraulic_permeability"))
    {
        DBUG("setting parameter `solid_hydraulic_permeability' to isotropic value `%g'", *par);
        const auto dim = BP::_mesh.getDimension();
        _assembly_params.solid_perm_tensor
                = Eigen::MatrixXd::Identity(dim, dim) * (*par);
    }

    // reactive system
    _assembly_params.react_sys = std::move(
        Ads::Adsorption::newInstance(config.getConfSubtree("reactive_system")));

    // matrix order
    {
        auto const order = config.getConfParam<std::string>("global_matrix_order");
        DBUG("global_matrix_order: %s", order.c_str());

        if (order == "BY_COMPONENT")
            _global_matrix_order = AssemblerLib::ComponentOrder::BY_COMPONENT;
        else if (order == "BY_LOCATION")
            _global_matrix_order = AssemblerLib::ComponentOrder::BY_LOCATION;
        else {
            ERR("unknown global matrix order `%s'", order.c_str());
            std::abort();
        }
    }

    // debug output
    if (auto const param = config.getConfParamOptional<bool>("output_element_matrices"))
    {
        DBUG("output_element_matrices: %s", (*param) ? "true" : "false");

        _assembly_params.output_element_matrices = *param;
    }

    // debug output
    if (auto const param = config.getConfParamOptional<bool>("output_iteration_results"))
    {
        DBUG("output_iteration_results: %s", (*param) ? "true" : "false");

        BP::_process_output.output_iteration_results = *param;
    }

    // debug output
    if (auto const param = config.getConfParamOptional<bool>("output_global_matrix"))
    {
        DBUG("output_global_matrix: %s", (*param) ? "true" : "false");

        BP::_process_output.output_global_matrix = *param;
    }
}

template<typename GlobalSetup>
void
TESProcess<GlobalSetup>::
createLocalAssemblers()
{
    switch (BP::_mesh.getDimension())
    {
    case 1: createLocalAssemblers<1>(); break;
    case 2: createLocalAssemblers<2>(); break;
    case 3: createLocalAssemblers<3>(); break;
    default:
        ERR("Invalid mesh dimension. Aborting.");
        std::abort();
    }
}

template<typename GlobalSetup>
template <unsigned GlobalDim>
void
TESProcess<GlobalSetup>::
createLocalAssemblers()
{
    DBUG("Create local assemblers.");
    // Populate the vector of local assemblers.
    _local_assemblers.resize(BP::_mesh.getNElements());
    // Shape matrices initializer
    using LocalDataInitializer = AssemblerLib::LocalDataInitializer<
        TES::TESLocalAssemblerInterface,
        TES::TESLocalAssembler,
        typename GlobalSetup::MatrixType,
        typename GlobalSetup::VectorType,
        GlobalDim,
        AssemblyParams const&>;

    LocalDataInitializer initializer;

    using LocalAssemblerBuilder =
        AssemblerLib::LocalAssemblerBuilder<
            MeshLib::Element,
            LocalDataInitializer>;

    LocalAssemblerBuilder local_asm_builder(
        initializer, *BP::_local_to_global_index_map);

    DBUG("Calling local assembler builder for all mesh elements.");
    BP::_global_setup.transform(
                local_asm_builder,
                BP::_mesh.getElements(),
                _local_assemblers,
                BP::_integration_order,
                _assembly_params);

    // TODO move somewhere else/make obsolete
    DBUG("Initialize TESProcess.");

    // for extrapolation of secondary variables
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> all_mesh_subsets_single_component;
    all_mesh_subsets_single_component.emplace_back(
                new MeshLib::MeshSubsets(BP::_mesh_subset_all_nodes));
    _local_to_global_index_map_single_component.reset(
                new AssemblerLib::LocalToGlobalIndexMap(
                    std::move(all_mesh_subsets_single_component), _global_matrix_order)
                );

    _extrapolator.reset(new ExtrapolatorImpl(*_local_to_global_index_map_single_component));
}

template<typename GlobalSetup>
void TESProcess<GlobalSetup>::
assembleConcreteProcess(
        const double t, GlobalVector const& x,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble TESProcess.");

    // Call global assembler for each local assembly item.
    BP::_global_setup.execute(*BP::_global_assembler,
                              _local_assemblers, t, x, M, K, b);

#ifndef NDEBUG
    if (_total_iteration == 0)
    {
        MathLib::BLAS::finalizeAssembly(M);
        MathLib::BLAS::finalizeAssembly(K);
        MathLib::BLAS::finalizeAssembly(b);

        // TODO [CL] Those files will be written to the working directory.
        //           Relative path needed.
        M.write("global_matrix_M.txt");
        K.write("global_matrix_K.txt");
        b.write("global_vector_b.txt");
    }
#endif
}

template<typename GlobalSetup>
void
TESProcess<GlobalSetup>::
preTimestep(GlobalVector const& /*x*/, const double t, const double delta_t)
{
    DBUG("new timestep");

    _assembly_params.delta_t = delta_t;
    _assembly_params.current_time = t;
    ++ _timestep; // TODO remove that
}

template<typename GlobalSetup>
void
TESProcess<GlobalSetup>::
preIteration(const unsigned iter, GlobalVector const& /*x*/)
{
    _assembly_params.iteration_in_current_timestep = iter;
}

template<typename GlobalSetup>
NumLib::IterationResult
TESProcess<GlobalSetup>::
postIteration(GlobalVector const& x)
{
    if (BP::_process_output.output_iteration_results)
    {
        DBUG("output results of iteration %li", _total_iteration);
        std::string fn = "tes_iter_" + std::to_string(_total_iteration) +
                         + "_ts_" + std::to_string(_timestep)
                         + "_" +    std::to_string(_assembly_params.iteration_in_current_timestep)
                         + ".vtu";

        BP::output(fn, 0, x);
    }

    bool check_passed = true;

    if (!Trafo::constrained)
    {
        // bounds checking only has to happen if the vapour mass fraction is non-logarithmic.

        auto do_check = [&](
                std::vector<double> const& local_x,
                AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& /*r_c_indices*/,
                LocalAssembler& loc_asm)
        {
            if (!loc_asm.checkBounds(local_x)) check_passed = false;
        };

        auto check_variable_bounds
        = [&](std::size_t id, LocalAssembler* const loc_asm)
        {
            BP::_global_assembler->passLocalVector(
                        do_check, id, x, *loc_asm);
        };

        // TODO Short-circuit evaluation that stops after the first error.
        //      But maybe that's not what I want to use here.
        BP::_global_setup.execute(
                    check_variable_bounds, _local_assemblers);
    }

    if (!check_passed)
        return NumLib::IterationResult::REPEAT_ITERATION;


    // TODO remove
    DBUG("ts %lu iteration %lu (%lu) try XXXXXX accepted", _timestep, _total_iteration,
         _assembly_params.iteration_in_current_timestep);

    ++ _assembly_params.iteration_in_current_timestep;
    ++_total_iteration;

    return NumLib::IterationResult::SUCCESS;
}

template<typename GlobalSetup>
typename TESProcess<GlobalSetup>::GlobalVector
TESProcess<GlobalSetup>::
computeVapourPartialPressure(typename TESProcess::GlobalVector const& x,
                             AssemblerLib::LocalToGlobalIndexMap const& dof_table)
{
    (void) x; (void) dof_table;
    return GlobalVector{};

#if 0
    case SecondaryVariables::VAPOUR_PARTIAL_PRESSURE:
    {
        IntegrationMethod_ integration_method(_integration_order);
        auto const n_integration_points = integration_method.getNPoints();

        auto& pVs = *_integration_point_values_cache;
        pVs.clear();
        pVs.reserve(n_integration_points);

        auto const& ps = nodal_dof.getElementNodalValues(0); // TODO [CL] use constants for DOF indices
        auto const& xs = nodal_dof.getElementNodalValues(2);

        auto const& AP = _data.getAssemblyParameters();

        for (auto const& sm : _shape_matrices)
        {
            double p, xm;

            using Array = std::array<double*, 1>;
            NumLib::shapeFunctionInterpolate(ps, sm.N, Array{ &p  });
            NumLib::shapeFunctionInterpolate(xs, sm.N, Array{ &xm });

            // TODO: Dalton's law method
            auto const xn = Ads::Adsorption::get_molar_fraction(xm, AP.M_react, AP.M_inert);
            pVs.push_back(p * xn);
        }

        return pVs;
    }
#endif
}

template<typename GlobalSetup>
typename TESProcess<GlobalSetup>::GlobalVector
TESProcess<GlobalSetup>::
computeRelativeHumidity(typename TESProcess::GlobalVector const& x,
                        AssemblerLib::LocalToGlobalIndexMap const& dof_table)
{
    (void) x; (void) dof_table;
    return GlobalVector{};

#if 0
    case SecondaryVariables::RELATIVE_HUMIDITY:
    {
        IntegrationMethod_ integration_method(_integration_order);
        auto const n_integration_points = integration_method.getNPoints();

        auto& rhs = *_integration_point_values_cache;
        rhs.clear();
        rhs.reserve(n_integration_points);

        auto const& nodal_vals = nodal_dof.getElementNodalValues();

        auto const& AP = _data.getAssemblyParameters();

        for (auto const& sm : _shape_matrices)
        {
            double p, T, xm;

            using Array = std::array<double*, 3>;
            NumLib::shapeFunctionInterpolate(nodal_vals, sm.N, Array{ &p, &T, &xm });

            // TODO: Dalton's law method
            auto const xn = Ads::Adsorption::get_molar_fraction(xm, AP.M_react, AP.M_inert);
            auto const pS = Ads::Adsorption::get_equilibrium_vapour_pressure(T);
            rhs.push_back(p * xn / pS);
        }

        return rhs;
    }
#endif
}

template<typename GlobalSetup>
typename TESProcess<GlobalSetup>::GlobalVector
TESProcess<GlobalSetup>::
computeEquilibriumLoading(typename TESProcess::GlobalVector const& x,
                          AssemblerLib::LocalToGlobalIndexMap const& dof_table)
{
    (void) x; (void) dof_table;
    return GlobalVector{};

#if 0
    case SecondaryVariables::EQUILIBRIUM_LOADING:
    {
        IntegrationMethod_ integration_method(_integration_order);
        auto const n_integration_points = integration_method.getNPoints();

        auto& Cs = *_integration_point_values_cache;
        Cs.clear();
        Cs.reserve(n_integration_points);

        auto const nodal_vals = nodal_dof.getElementNodalValues();

        auto const& AP = _data.getAssemblyParameters();

        for (auto const& sm : _shape_matrices)
        {
            double p, T, xm;

            using Array = std::array<double*, 3>;
            NumLib::shapeFunctionInterpolate(nodal_vals, sm.N, Array{ &p, &T, &xm });

            // TODO: Dalton's law method
            auto const xn = Ads::Adsorption::get_molar_fraction(xm, AP.M_react, AP.M_inert);
            auto const pV = p * xn;
            if (pV < 0.0) {
                Cs.push_back(0.0);
            } else {
                Cs.push_back(AP.react_sys->get_equilibrium_loading(pV, T, AP.M_react));
            }
        }

        return Cs;
    }
#endif
}

template<typename GlobalSetup>
SecondaryVariableFunctions<typename TESProcess<GlobalSetup>::GlobalVector>
TESProcess<GlobalSetup>::
makeExtrapolator(TESIntPtVariables const var) const
{
    auto const eval_field = [var, this](
            GlobalVector const& x,
            AssemblerLib::LocalToGlobalIndexMap const& /*dof_table*/
            ) -> GlobalVector
    {
        _extrapolator->extrapolate(x, _local_assemblers, var);
        return _extrapolator->getNodalValues();
    };

    auto const eval_residuals = [var, this](
            GlobalVector const& x,
            AssemblerLib::LocalToGlobalIndexMap const& /*dof_table*/
            ) -> GlobalVector
    {
        _extrapolator->calculateResiduals(x, _local_assemblers, var);
        return _extrapolator->getElementResiduals();
    };
    return { eval_field, eval_residuals };
}

template<typename GlobalSetup>
TESProcess<GlobalSetup>::
~TESProcess()
{
    for (auto p : _local_assemblers)
        delete p;
}


// Explicitly instantiate TESProcess for GlobalSetupType.
template class TESProcess<GlobalSetupType>;

} // namespace TES

}   // namespace ProcessLib
