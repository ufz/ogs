/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_TESPROCESS_IMPL_H_
#define PROCESS_LIB_TESPROCESS_IMPL_H_

#include <cassert>
#include <cstdio>

#include <logog/include/logog.hpp>

#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "NumLib/TimeStepping/Algorithms/FixedTimeStepping.h"
#include "MathLib/LinAlg/VectorNorms.h"

#include "TESProcess.h"


namespace
{

template <class ConfigTree>
static ProcessLib::ProcessVariable const*
find_variable(ConfigTree const& config,
              std::string const& variable_role,
              std::vector<ProcessLib::ProcessVariable> const& variables)
{
    std::string const name = config.template getConfParam<std::string>(variable_role);

    auto variable = std::find_if(variables.cbegin(), variables.cend(),
            [&name](ProcessLib::ProcessVariable const& v) {
                    DBUG("proc var `%s'", v.getName().c_str());
                return v.getName() == name;
            });

    if (variable == variables.end())
    {
        ERR("Expected process variable \'%s\' not found in provided variables list.",
            name.c_str());
        std::abort();
    }


    DBUG("Found process variable %s for role %s.",
         name.c_str(), variable_role.c_str());

    return &*variable;
}

} // anonymous namespace


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
        const std::string vars[NODAL_DOF] = { "fluid_pressure",
                                              "temperature",
                                              "vapour_mass_fraction" };

        auto const proc_vars = config.getConfSubtree("process_variables");

        for (unsigned i=0; i<NODAL_DOF; ++i)
        {
            auto variable = find_variable(proc_vars, vars[i], variables);

            BP::_process_variables[i] = *const_cast<ProcessVariable*>(variable);
        }
    }

    // secondary variables
    if (auto proc_vars = config.getConfSubtreeOptional("secondary_variables"))
    {
        auto add_secondary_variable =
                [this, &proc_vars](
                std::string const& var, SecondaryVariables type, unsigned num_components)
                -> void
        {
            if (auto variable = proc_vars->getConfParamOptional<std::string>(var))
            {
                _secondary_process_vars.emplace_back(type, *variable, num_components);
            }
        };

        add_secondary_variable("solid_density", SecondaryVariables::SOLID_DENSITY, 1);
        add_secondary_variable("reaction_rate", SecondaryVariables::REACTION_RATE, 1);
        add_secondary_variable("velocity_x",    SecondaryVariables::VELOCITY_X,    1);
        if (BP::_mesh.getDimension() >= 2) add_secondary_variable("velocity_y",    SecondaryVariables::VELOCITY_Y,    1);
        if (BP::_mesh.getDimension() >= 3) add_secondary_variable("velocity_z",    SecondaryVariables::VELOCITY_Z,    1);

        add_secondary_variable("vapour_partial_pressure", SecondaryVariables::VAPOUR_PARTIAL_PRESSURE, 1);
        add_secondary_variable("relative_humidity",       SecondaryVariables::RELATIVE_HUMIDITY,       1);
        add_secondary_variable("loading",                 SecondaryVariables::LOADING,                 1);
        add_secondary_variable("equilibrium_loading",     SecondaryVariables::EQUILIBRIUM_LOADING,     1);
        add_secondary_variable("reaction_damping_factor", SecondaryVariables::REACTION_DAMPING_FACTOR, 1);
    }

    // variables for output
    if (auto output = config.getConfSubtreeOptional("output")) {
        if (auto out_vars = output->getConfSubtreeOptional("variables"))
        {
            for (auto out_var : out_vars->getConfParamList<std::string>("variable"))
            {
                if (_output_variables.find(out_var) != _output_variables.cend())
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
                    auto pred2 = [&out_var](std::tuple<SecondaryVariables, std::string, unsigned> const& p) {
                        return std::get<1>(p) == out_var;
                    };

                    // check if secondary variable
                    auto const& pcs_var2 = std::find_if(
                        _secondary_process_vars.cbegin(), _secondary_process_vars.cend(),
                        pred2);

                    if (pcs_var2 == _secondary_process_vars.cend())
                    {
                        ERR("Output variable `%s' is neither a process variable nor a"
                            " secondary variable", out_var.c_str());
                        std::abort();
                    }
                }

                DBUG("adding output variable `%s'", out_var.c_str());
                _output_variables.insert(out_var);
            }

            if (auto out_resid = output->getConfParamOptional<bool>("output_extrapolation_residuals"))
            {
                _output_residuals = *out_resid;
            }
        }
    }

    {
        std::vector<std::pair<const std::string, double*> > params{
            { "fluid_specific_heat_source",            &_assembly_params._fluid_specific_heat_source },
            { "fluid_specific_isobaric_heat_capacity", &_assembly_params._cpG },
            { "solid_specific_heat_source",            &_assembly_params._solid_specific_heat_source },
            { "solid_heat_conductivity",               &_assembly_params._solid_heat_cond },
            { "solid_specific_isobaric_heat_capacity", &_assembly_params._cpS },
            { "tortuosity",                            &_assembly_params._tortuosity },
            { "diffusion_coefficient",                 &_assembly_params._diffusion_coefficient_component },
            { "porosity",                              &_assembly_params._poro },
            { "solid_density_dry",                     &_assembly_params._rho_SR_dry },
            { "solid_density_initial",                 &_assembly_params._initial_solid_density }
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
        _assembly_params._solid_perm_tensor
                = Eigen::MatrixXd::Identity(dim, dim) * (*par);
    }

    // reactive system
    _assembly_params._reaction_system = std::move(
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

        _assembly_params._output_element_matrices = *param;
    }

    // debug output
    if (auto const param = config.getConfParamOptional<bool>("output_iteration_results"))
    {
        DBUG("output_iteration_results: %s", (*param) ? "true" : "false");

        _output_iteration_results = *param;
    }

    // debug output
    if (auto const param = config.getConfParamOptional<bool>("output_global_matrix"))
    {
        DBUG("output_global_matrix: %s", (*param) ? "true" : "false");

        _output_global_matrix = *param;
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
        TES::LocalAssemblerDataInterface,
        TES::LocalAssemblerData,
        typename GlobalSetup::MatrixType,
        typename GlobalSetup::VectorType,
        GlobalDim>;

    LocalDataInitializer initializer;

    using LocalAssemblerBuilder =
        AssemblerLib::LocalAssemblerBuilder<
            MeshLib::Element,
            LocalDataInitializer>;

    LocalAssemblerBuilder local_asm_builder(
        initializer, *BP::_local_to_global_index_map);

    DBUG("Calling local assembler builder for all mesh elements.");
    _global_setup.transform(
                local_asm_builder,
                BP::_mesh.getElements(),
                _local_assemblers,
                _integration_order,
                this);

    DBUG("Create global assembler.");
    _global_assembler.reset(
        new GlobalAssembler(*BP::_local_to_global_index_map));


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
    const double current_time = 0.0;
    DBUG("Assemble TESProcess.");

    // _assembly_params._delta_t = delta_t; // TODO fix
    _assembly_params._iteration_in_current_timestep = 0;
    _assembly_params._current_time = current_time;
    ++ _timestep;


    // from singlePicardIteration()

    bool iteration_accepted = false;
    unsigned num_try = 0;

    do
    {
        INFO("-> TES process try number %u in current picard iteration", num_try);
        _assembly_params._number_of_try_of_iteration = num_try;

        // Call global assembler for each local assembly item.
        _global_setup.execute(*BP::_global_assembler, _local_assemblers,
                              t, x, M, K, b);

#if 0 && defined(OGS_USE_EIGEN) && ! defined(OGS_USE_EIGENLIS)
        // TODO put that somewhere
        MathLib::scaleDiagonal(*_A, *_rhs);
#endif

#if 0 && defined(OGS_USE_EIGENLIS)
        // TODO put that somewhere

        // scaling
        typename GlobalMatrix::RawMatrixType AT = _A->getRawMatrix().transpose();

        for (unsigned dof = 0; dof < NODAL_DOF; ++dof)
        {
            auto const& trafo = (dof == 0) ? _assembly_params.trafo_p
                              : (dof == 1) ? _assembly_params.trafo_T
                                           : _assembly_params.trafo_x;

            for (std::size_t i = 0; i < BP::_mesh.getNNodes(); ++i)
            {
                MeshLib::Location loc(BP::_mesh.getID(), MeshLib::MeshItemType::Node, i);
                auto const idx = BP::_local_to_global_index_map->getGlobalIndex(loc, dof);

                AT.row(idx) *= trafo.dxdy(0);
                x_curr[idx] /= trafo.dxdy(0);
            }
        }

        _A->getRawMatrix() = AT.transpose();
#endif

#ifndef NDEBUG
        if (_total_iteration == 0 && num_try == 0 && _output_global_matrix)
        {
            // TODO fix after PETSc
            // MathLib::BLAS::finalizeAssembly(*_A); //  MathLib::finalizeMatrixAssembly(*_A);

            // TODO fix
#if 0
            // TODO [CL] Those files will be written to the working directory.
            //           Relative path needed.
            _A->write("global_matrix.txt");
            _rhs->write("global_rhs.txt");
#endif
        }
#endif

#if 0 && defined(OGS_USE_EIGENLIS)
        // TODO put that somewhere

        // scale back
        for (unsigned dof = 0; dof < NODAL_DOF; ++dof)
        {
            auto const& trafo = (dof == 0) ? _assembly_params.trafo_p
                              : (dof == 1) ? _assembly_params.trafo_T
                                           : _assembly_params.trafo_x;

            for (std::size_t i = 0; i < BP::_mesh.getNNodes(); ++i)
            {
                MeshLib::Location loc(BP::_mesh.getID(), MeshLib::MeshItemType::Node, i);
                auto const idx = BP::_local_to_global_index_map->getGlobalIndex(loc, dof);

                // TODO: _A
                x_curr[idx] *= trafo.dxdy(0);
            }
        }
#endif

        if (_output_iteration_results)
        {
            DBUG("output results of iteration %li", _total_iteration);
            std::string fn = "tes_iter_" + std::to_string(_total_iteration) +
                             + "_ts_" + std::to_string(_timestep)
                             + "_" +    std::to_string(_assembly_params._iteration_in_current_timestep)
                             + "_" +    std::to_string(num_try)
                             + ".vtu";

            postTimestep(fn, 0);
        }

        bool check_passed = true;

        // TODO put to post timestep.
        if (!Trafo::constrained)
        {
            // bounds checking only has to happen if the vapour mass fraction is non-logarithmic.

            auto& ga = *_global_assembler;

            auto check_variable_bounds
            = [&ga, &check_passed](
              std::size_t id, LocalAssembler* const loc_asm)
            {
                // DBUG("%lu", id);

                std::vector<double> const* localX;
                std::vector<double> const* localX_pts;

                // TODO fix
                // ga.getLocalNodalValues(id, localX, localX_pts);

                // if (!loc_asm->checkBounds(*localX, *localX_pts)) check_passed = false;
            };

            _global_setup.execute(check_variable_bounds, _local_assemblers);

            if (!check_passed)
            {
                // TODO fix
                // x_curr = x_prev_iter;
            }
        }

        iteration_accepted = check_passed;

        ++num_try;
    }
    while(! iteration_accepted);

    DBUG("ts %lu iteration %lu (%lu) try %u accepted", _timestep, _total_iteration,
         _assembly_params._iteration_in_current_timestep, num_try-1);

    ++ _assembly_params._iteration_in_current_timestep;
    ++_total_iteration;
}


template<typename GlobalSetup>
void
TESProcess<GlobalSetup>::
postTimestep(const std::string& file_name, const unsigned /*timestep*/)
// TODO [CL] remove second parameter
{
    INFO("postprocessing timestep");

    /*
    std::puts("---- solution ----");
    printGlobalVector(_x->getRawVector());
    */

    auto count = [](MeshLib::Mesh const& mesh, MeshLib::MeshItemType type) -> std::size_t
    {
        switch (type) {
        case MeshLib::MeshItemType::Cell: return mesh.getNElements();
        case MeshLib::MeshItemType::Node: return mesh.getNNodes();
        default: break;
        }
        return 0;
    };

    auto get_or_create_mesh_property = [this, &count](std::string const& property_name, MeshLib::MeshItemType type)
    {
        // Get or create a property vector for results.
        boost::optional<MeshLib::PropertyVector<double>&> result;

        auto const N = count(BP::_mesh, type);

        if (BP::_mesh.getProperties().hasPropertyVector(property_name))
        {
            result = BP::_mesh.getProperties().template
                getPropertyVector<double>(property_name);
        }
        else
        {
            result = BP::_mesh.getProperties().template
                createNewPropertyVector<double>(property_name, type);
            result->resize(N);
        }
        assert(result && result->size() == N);

        return result;
    };

    auto add_primary_var = [this, &get_or_create_mesh_property](const unsigned vi)
    {
        std::string const& property_name = BP::_process_variables[vi].get().getName();
        if (_output_variables.find(property_name) == _output_variables.cend())
            return;

        DBUG("  process var %s", property_name.c_str());

        auto result = get_or_create_mesh_property(property_name, MeshLib::MeshItemType::Node);
        assert(result->size() == BP::_mesh.getNNodes());

        // Copy result
        for (std::size_t i = 0; i < BP::_mesh.getNNodes(); ++i)
        {
            MeshLib::Location loc(BP::_mesh.getID(), MeshLib::MeshItemType::Node, i);
            auto const idx = BP::_local_to_global_index_map->getGlobalIndex(loc, vi);
            assert(!std::isnan((*_x)[idx]));
            (*result)[i] = (*_x)[idx];
        }
    };

    assert(_x->size() == NODAL_DOF * BP::_mesh.getNNodes());
    for (unsigned vi=0; vi!=NODAL_DOF; ++vi)
    {
        add_primary_var(vi);
    }


    auto add_secondary_var = [this, &get_or_create_mesh_property]
                             (SecondaryVariables const property,
                             std::string const& property_name,
                             const unsigned num_components
                             ) -> void
    {
        assert(num_components == 1); // TODO [CL] implement other cases
        (void) num_components;

        {
            if (_output_variables.find(property_name) == _output_variables.cend())
                return;

            DBUG("  process var %s", property_name.c_str());

            auto result = get_or_create_mesh_property(property_name, MeshLib::MeshItemType::Node);
            assert(result->size() == BP::_mesh.getNNodes());

            _extrapolator->extrapolate(*_x, *BP::_local_to_global_index_map, _local_assemblers, property);
            auto const& nodal_values = _extrapolator->getNodalValues();

            // Copy result
            for (std::size_t i = 0; i < BP::_mesh.getNNodes(); ++i)
            {
                assert(!std::isnan(nodal_values[i]));
                (*result)[i] = nodal_values[i];
            }
        }

        if (_output_residuals) {
            DBUG("  process var %s residual", property_name.c_str());
            auto const& property_name_res = property_name + "_residual";

            auto result = get_or_create_mesh_property(property_name_res, MeshLib::MeshItemType::Cell);
            assert(result->size() == BP::_mesh.getNElements());

            _extrapolator->calculateResiduals(*_x, *BP::_local_to_global_index_map, _local_assemblers, property);
            auto const& residuals = _extrapolator->getElementResiduals();

            // Copy result
            for (std::size_t i = 0; i < BP::_mesh.getNElements(); ++i)
            {
                assert(!std::isnan(residuals[i]));
                (*result)[i] = residuals[i];
            }
        }
    };

    for (auto const& p : _secondary_process_vars)
    {
        add_secondary_var(std::get<0>(p), std::get<1>(p), std::get<2>(p));
    }


    // Write output file
    FileIO::VtuInterface vtu_interface(&this->_mesh, vtkXMLWriter::Binary, true);
    vtu_interface.writeToFile(file_name);
}


template<typename GlobalSetup>
void
TESProcess<GlobalSetup>::
post(std::string const& file_name)
{
    DBUG("Postprocessing TESProcess.");
    std::string const property_name = "Result";

    // Get or create a property vector for results.
    boost::optional<MeshLib::PropertyVector<double>&> result;
    if (BP::_mesh.getProperties().hasPropertyVector(property_name))
    {
        result = BP::_mesh.getProperties().template
            getPropertyVector<double>(property_name);
    }
    else
    {
        result = BP::_mesh.getProperties().template
            createNewPropertyVector<double>(property_name,
                MeshLib::MeshItemType::Node);
        result->resize(_x->size());
    }
    assert(result && result->size() == _x->size());

    // Copy result
    for (std::size_t i = 0; i < _x->size(); ++i)
        (*result)[i] = (*_x)[i];

    // Write output file
    FileIO::VtuInterface vtu_interface(&this->_mesh, vtkXMLWriter::Binary, true);
    vtu_interface.writeToFile(file_name);
}

template<typename GlobalSetup>
TESProcess<GlobalSetup>::
~TESProcess()
{
    for (auto p : _local_assemblers)
        delete p;
}

} // namespace TES

} // namespace ProcessLib

#endif  // PROCESS_LIB_TESPROCESS_IMPL_H_
