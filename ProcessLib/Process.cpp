/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Process.h"

#include "BaseLib/Functional.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "NumLib/ODESolver/ConvergenceCriterionPerComponent.h"
#include "GlobalVectorFromNamedFunction.h"
#include "ProcessVariable.h"

namespace ProcessLib
{
Process::Process(
    MeshLib::Mesh& mesh,
    NonlinearSolver& nonlinear_solver,
    std::unique_ptr<TimeDiscretization>&& time_discretization,
    std::unique_ptr<NumLib::ConvergenceCriterion>&& convergence_criterion,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    SecondaryVariableCollection&& secondary_variables,
    ProcessOutput&& process_output,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : _mesh(mesh),
      _secondary_variables(std::move(secondary_variables)),
      _process_output(std::move(process_output)),
      _named_function_caller(std::move(named_function_caller)),
      _nonlinear_solver(nonlinear_solver),
      _time_discretization(std::move(time_discretization)),
      _convergence_criterion(std::move(convergence_criterion)),
      _process_variables(std::move(process_variables)),
      _boundary_conditions(parameters)
{
}

void Process::output(std::string const& file_name,
                     const unsigned /*timestep*/,
                     GlobalVector const& x) const
{
    doProcessOutput(file_name, x, _mesh, *_local_to_global_index_map,
                    _process_variables, _secondary_variables, _process_output);
}

void Process::initialize()
{
    DBUG("Initialize process.");

    DBUG("Construct dof mappings.");
    constructDofTable();

    DBUG("Compute sparsity pattern");
    computeSparsityPattern();

    DBUG("Initialize the extrapolator");
    initializeExtrapolator();

    initializeConcreteProcess(*_local_to_global_index_map, _mesh,
                              _integration_order);

    finishNamedFunctionsInitialization();

    DBUG("Initialize boundary conditions.");
    _boundary_conditions.addBCsForProcessVariables(
        _process_variables, *_local_to_global_index_map, _integration_order);
}

void Process::setInitialConditions(double const t, GlobalVector& x)
{
    DBUG("Set initial conditions.");
    std::size_t const n_nodes = _mesh.getNumberOfNodes();

    SpatialPosition pos;

    for (int variable_id = 0;
         variable_id < static_cast<int>(_process_variables.size());
         ++variable_id)
    {
        ProcessVariable& pv = _process_variables[variable_id];
        auto const& ic = pv.getInitialCondition();

        auto const num_comp = pv.getNumberOfComponents();

        for (std::size_t node_id = 0; node_id < n_nodes; ++node_id)
        {
            MeshLib::Location const l(_mesh.getID(),
                                      MeshLib::MeshItemType::Node, node_id);

            pos.setNodeID(node_id);
            auto const& tup = ic.getTuple(t, pos);

            for (int comp_id = 0; comp_id < num_comp; ++comp_id)
            {
                auto global_index =
                    std::abs(_local_to_global_index_map->getGlobalIndex(
                        l, variable_id, comp_id));
#ifdef USE_PETSC
                // The global indices of the ghost entries of the global matrix
                // or the global vectors need to be set as negative values for
                // equation assembly, however the global indices start from
                // zero. Therefore, any ghost entry with zero index is assigned
                // an negative value of the vector size or the matrix dimension.
                // To assign the initial value for the ghost entries, the
                // negative indices of the ghost entries are restored to zero.
                if (global_index == x.size())
                    global_index = 0;
#endif
                x.set(global_index, tup[comp_id]);
            }
        }
    }
}

MathLib::MatrixSpecifications Process::getMatrixSpecifications() const
{
    auto const& l = *_local_to_global_index_map;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern};
}

void Process::assemble(const double t, GlobalVector const& x, GlobalMatrix& M,
                       GlobalMatrix& K, GlobalVector& b)
{
    assembleConcreteProcess(t, x, M, K, b);
    _boundary_conditions.apply(t, x, K, b);
}

void Process::assembleJacobian(const double t, GlobalVector const& x,
                               GlobalVector const& xdot, const double dxdot_dx,
                               GlobalMatrix const& M, const double dx_dx,
                               GlobalMatrix const& K, GlobalMatrix& Jac)
{
    assembleJacobianConcreteProcess(t, x, xdot, dxdot_dx, M, dx_dx, K, Jac);

    // TODO In this method one could check if the user wants to use an
    //      analytical or a numerical Jacobian. Then the right
    //      assembleJacobianConcreteProcess() method will be chosen.
    //      Additionally in the default implementation of said method one
    //      could provide a fallback to a numerical Jacobian. However, that
    //      would be in a sense implicit behaviour and it might be better to
    //      just abort, as is currently the case.
    //      In order to implement the Jacobian assembly entirely, in addition
    //      to the operator() in VectorMatrixAssembler there has to be a method
    //      that dispatches the Jacobian assembly.
    //      Similarly, then the NeumannBC specialization of VectorMatrixAssembler
    //      probably can be merged into the main class s.t. one has only one
    //      type of VectorMatrixAssembler (for each equation type) with the
    //      three methods assemble(), assembleJacobian() and assembleNeumannBC().
    //      That list can be extended, e.g. by methods for the assembly of
    //      source terms.
    //      UPDATE: Probably it is better to keep a separate NeumannBC version of the
    //      VectoMatrixAssembler since that will work for all kinds of processes.
}

void Process::assembleJacobianConcreteProcess(
    const double /*t*/, GlobalVector const& /*x*/, GlobalVector const& /*xdot*/,
    const double /*dxdot_dx*/, GlobalMatrix const& /*M*/,
    const double /*dx_dx*/, GlobalMatrix const& /*K*/, GlobalMatrix& /*Jac*/)
{
    OGS_FATAL(
        "The concrete implementation of this Process did not override the"
        " assembleJacobianConcreteProcess() method."
        " Hence, no analytical Jacobian is provided for this process"
        " and the Newton-Raphson method cannot be used to solve it.");
}

void Process::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes.reset(
        new MeshLib::MeshSubset(_mesh, &_mesh.getNodes()));

    // Collect the mesh subsets in a vector.
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> all_mesh_subsets;
    for (ProcessVariable const& pv : _process_variables)
    {
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            pv.getNumberOfComponents(),
            [&]() {
                return std::unique_ptr<MeshLib::MeshSubsets>{
                    new MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()}};
            });
    }

    _local_to_global_index_map.reset(new NumLib::LocalToGlobalIndexMap(
        std::move(all_mesh_subsets), NumLib::ComponentOrder::BY_LOCATION));

    if (auto* conv_crit =
            dynamic_cast<NumLib::ConvergenceCriterionPerComponent*>(
                _convergence_criterion.get())) {
        conv_crit->setDOFTable(*_local_to_global_index_map, _mesh);
    }
}

void Process::initializeExtrapolator()
{
    NumLib::LocalToGlobalIndexMap const* dof_table_single_component;
    bool manage_storage;

    if (_local_to_global_index_map->getNumberOfComponents() == 1)
    {
        // For single-variable-single-component processes reuse the existing DOF
        // table.
        dof_table_single_component = _local_to_global_index_map.get();
        manage_storage = false;
    }
    else
    {
        // Otherwise construct a new DOF table.
        std::vector<std::unique_ptr<MeshLib::MeshSubsets>>
            all_mesh_subsets_single_component;
        all_mesh_subsets_single_component.emplace_back(
            new MeshLib::MeshSubsets(_mesh_subset_all_nodes.get()));

        dof_table_single_component = new NumLib::LocalToGlobalIndexMap(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);
        manage_storage = true;
    }

    std::unique_ptr<NumLib::Extrapolator> extrapolator(
        new NumLib::LocalLinearLeastSquaresExtrapolator(
            *dof_table_single_component));

    // TODO Later on the DOF table can change during the simulation!
    _extrapolator_data = ExtrapolatorData(
        std::move(extrapolator), dof_table_single_component, manage_storage);
}

void Process::finishNamedFunctionsInitialization()
{
    _named_function_caller.applyPlugs();

    for (auto const& named_function :
         _named_function_caller.getNamedFunctions()) {
        auto const& name = named_function.getName();
        // secondary variables generated from named functions have the prefix
        // "fct_".
        _secondary_variables.addSecondaryVariable(
            "fct_" + name, 1,
            {BaseLib::easyBind(
                 &GlobalVectorFromNamedFunction::call,
                 GlobalVectorFromNamedFunction(
                     _named_function_caller.getSpecificFunctionCaller(name), _mesh,
                     getSingleComponentDOFTable(),
                     _secondary_variable_context)),
             nullptr});
    }
}

void Process::computeSparsityPattern()
{
    _sparsity_pattern =
        NumLib::computeSparsityPattern(*_local_to_global_index_map, _mesh);
}

void Process::preIteration(const unsigned iter, const GlobalVector &x)
{
    // In every new iteration cached values of secondary variables are expired.
    for (auto& cached_var : _cached_secondary_variables) {
        cached_var->expire();
    }

    preIterationConcreteProcess(iter, x);
}

NumLib::IterationResult Process::postIteration(const GlobalVector &x)
{
    return postIterationConcreteProcess(x);
}

}  // namespace ProcessLib
