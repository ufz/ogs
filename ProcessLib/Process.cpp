/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Process.h"

#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "NumLib/ODESolver/ConvergenceCriterionPerComponent.h"
#include "ParameterLib/Parameter.h"

#include "ProcessVariable.h"
#include "CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
Process::Process(
    std::string name_,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    SecondaryVariableCollection&& secondary_variables,
    const bool use_monolithic_scheme)
    : name(std::move(name_)),
      mesh_(mesh),
      secondary_variables_(std::move(secondary_variables)),
      global_assembler_(std::move(jacobian_assembler)),
      use_monolithic_scheme_(use_monolithic_scheme),
      coupled_solutions_(nullptr),
      integration_order_(integration_order),
      process_variables_(std::move(process_variables)),
      boundary_conditions_([&](const std::size_t number_of_process_variables)
                               -> std::vector<BoundaryConditionCollection> {
          std::vector<BoundaryConditionCollection> pcs_BCs;
          pcs_BCs.reserve(number_of_process_variables);
          for (std::size_t i = 0; i < number_of_process_variables; i++)
          {
              pcs_BCs.emplace_back(BoundaryConditionCollection(parameters));
          }
          return pcs_BCs;
      }(process_variables_.size())),
      source_term_collections_([&](const std::size_t number_of_processes)
                                   -> std::vector<SourceTermCollection> {
          std::vector<SourceTermCollection> pcs_sts;
          pcs_sts.reserve(number_of_processes);
          for (std::size_t i = 0; i < number_of_processes; i++)
          {
              pcs_sts.emplace_back(SourceTermCollection(parameters));
          }
          return pcs_sts;
      }(process_variables_.size()))
{
}

void Process::initializeProcessBoundaryConditionsAndSourceTerms(
    const NumLib::LocalToGlobalIndexMap& dof_table, const int process_id)
{
    auto const& per_process_variables = process_variables_[process_id];
    auto& per_process_BCs = boundary_conditions_[process_id];

    per_process_BCs.addBCsForProcessVariables(per_process_variables, dof_table,
                                              integration_order_, *this);

    auto& per_process_sts = source_term_collections_[process_id];
    per_process_sts.addSourceTermsForProcessVariables(
        per_process_variables, dof_table, integration_order_);
}

void Process::initializeBoundaryConditions()
{
    // The number of processes is identical to the size of process_variables_,
    // the vector contains variables for different processes. See the
    // documentation of process_variables_.
    const std::size_t number_of_processes = process_variables_.size();
    for (std::size_t pcs_id = 0; pcs_id < number_of_processes; pcs_id++)
    {
        initializeProcessBoundaryConditionsAndSourceTerms(
            *local_to_global_index_map_, pcs_id);
    }
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

    initializeConcreteProcess(*local_to_global_index_map_, mesh_,
                              integration_order_);

    DBUG("Initialize boundary conditions.");
    initializeBoundaryConditions();
}

void Process::setInitialConditions(const int process_id, double const t,
                                   GlobalVector& x)
{
    // getDOFTableOfProcess can be overloaded by the specific process.
    auto const& dof_table_of_process = getDOFTable(process_id);

    auto const& per_process_variables = process_variables_[process_id];
    for (std::size_t variable_id = 0;
         variable_id < per_process_variables.size();
         variable_id++)
    {
        MathLib::LinAlg::setLocalAccessibleVector(x);
        ParameterLib::SpatialPosition pos;

        auto const& pv = per_process_variables[variable_id];
        DBUG("Set the initial condition of variable {:s} of process {:d}.",
             pv.get().getName().data(), process_id);

        auto const& ic = pv.get().getInitialCondition();

        auto const num_comp = pv.get().getNumberOfComponents();

        for (int component_id = 0; component_id < num_comp; ++component_id)
        {
            auto const& mesh_subset =
                dof_table_of_process.getMeshSubset(variable_id, component_id);
            auto const mesh_id = mesh_subset.getMeshID();
            for (auto const* node : mesh_subset.getNodes())
            {
                MeshLib::Location const l(mesh_id, MeshLib::MeshItemType::Node,
                                          node->getID());

                pos.setNodeID(node->getID());
                auto const& ic_value = ic(t, pos);

                auto global_index =
                    std::abs(dof_table_of_process.getGlobalIndex(l, variable_id,
                                                                 component_id));
#ifdef USE_PETSC
                // The global indices of the ghost entries of the global
                // matrix or the global vectors need to be set as negative
                // values for equation assembly, however the global indices
                // start from zero. Therefore, any ghost entry with zero
                // index is assigned an negative value of the vector size
                // or the matrix dimension. To assign the initial value for
                // the ghost entries, the negative indices of the ghost
                // entries are restored to zero.
                if (global_index == x.size())
                    global_index = 0;
#endif
                x.set(global_index, ic_value[component_id]);
            }
        }
    }
    setInitialConditionsConcreteProcess(x, t);
}

MathLib::MatrixSpecifications Process::getMatrixSpecifications(
    const int /*process_id*/) const
{
    auto const& l = *local_to_global_index_map_;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &sparsity_pattern_};
}

void Process::updateDeactivatedSubdomains(double const time,
                                          const int process_id)
{
    auto const& variables_per_process = getProcessVariables(process_id);
    for (auto const& variable : variables_per_process)
    {
        variable.get().updateDeactivatedSubdomains(time);
    }
}

void Process::preAssemble(const double t, double const dt,
                          GlobalVector const& x)
{
    preAssembleConcreteProcess(t, dt, x);
}

void Process::assemble(const double t, double const dt,
                       std::vector<GlobalVector*> const& x,
                       std::vector<GlobalVector*> const& xdot,
                       int const process_id, GlobalMatrix& M, GlobalMatrix& K,
                       GlobalVector& b)
{
    MathLib::LinAlg::setLocalAccessibleVector(*x[process_id]);
    MathLib::LinAlg::setLocalAccessibleVector(*xdot[process_id]);

    assembleConcreteProcess(t, dt, x, xdot, process_id, M, K, b);

    // the last argument is for the jacobian, nullptr is for a unused jacobian
    boundary_conditions_[process_id].applyNaturalBC(t, x, process_id, K, b,
                                                    nullptr);

    // the last argument is for the jacobian, nullptr is for a unused jacobian
    source_term_collections_[process_id].integrate(t, *x[process_id], b,
                                                   nullptr);
}

void Process::assembleWithJacobian(const double t, double const dt,
                                   std::vector<GlobalVector*> const& x,
                                   GlobalVector const& xdot,
                                   const double dxdot_dx, const double dx_dx,
                                   int const process_id, GlobalMatrix& M,
                                   GlobalMatrix& K, GlobalVector& b,
                                   GlobalMatrix& Jac)
{
    MathLib::LinAlg::setLocalAccessibleVector(*x[process_id]);
    MathLib::LinAlg::setLocalAccessibleVector(xdot);

    assembleWithJacobianConcreteProcess(t, dt, x, xdot, dxdot_dx, dx_dx,
                                        process_id, M, K, b, Jac);

    // TODO: apply BCs to Jacobian.
    boundary_conditions_[process_id].applyNaturalBC(t, x, process_id, K, b,
                                                    &Jac);

    // the last argument is for the jacobian, nullptr is for a unused jacobian
    source_term_collections_[process_id].integrate(t, *x[process_id], b, &Jac);
}

void Process::constructDofTable()
{
    if (use_monolithic_scheme_)
    {
        constructMonolithicProcessDofTable();

        return;
    }

    // For staggered scheme:
    const int specified_prosess_id = 0;
    constructDofTableOfSpecifiedProsessStaggerdScheme(specified_prosess_id);
}

void Process::constructMonolithicProcessDofTable()
{
    // Create single component dof in every of the mesh nodes.
    mesh_subset_all_nodes_ =
        std::make_unique<MeshLib::MeshSubset>(mesh_, mesh_.getNodes());

    // Vector of mesh subsets.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;

    // Vector of the number of variable components
    std::vector<int> vec_var_n_components;
    // Collect the mesh subsets in a vector for the components of each
    // variables.
    for (ProcessVariable const& pv : process_variables_[0])
    {
        std::generate_n(std::back_inserter(all_mesh_subsets),
                        pv.getNumberOfComponents(),
                        [&]() { return *mesh_subset_all_nodes_; });
    }

    // Create a vector of the number of variable components
    for (ProcessVariable const& pv : process_variables_[0])
    {
        vec_var_n_components.push_back(pv.getNumberOfComponents());
    }

    local_to_global_index_map_ =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), vec_var_n_components,
            NumLib::ComponentOrder::BY_LOCATION);

    assert(local_to_global_index_map_);
}

void Process::constructDofTableOfSpecifiedProsessStaggerdScheme(
    const int specified_prosess_id)
{
    // Create single component dof in every of the mesh nodes.
    mesh_subset_all_nodes_ =
        std::make_unique<MeshLib::MeshSubset>(mesh_, mesh_.getNodes());

    // Vector of mesh subsets.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets;

    // Vector of the number of variable components
    std::vector<int> vec_var_n_components;
    // Collect the mesh subsets in a vector for each variables' components.
    std::generate_n(std::back_inserter(all_mesh_subsets),
                    process_variables_[specified_prosess_id][0]
                        .get()
                        .getNumberOfComponents(),
                    [&]() { return *mesh_subset_all_nodes_; });

    // Create a vector of the number of variable components.
    vec_var_n_components.push_back(process_variables_[specified_prosess_id][0]
                                       .get()
                                       .getNumberOfComponents());
    local_to_global_index_map_ =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), vec_var_n_components,
            NumLib::ComponentOrder::BY_LOCATION);

    assert(local_to_global_index_map_);
}

std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
Process::getDOFTableForExtrapolatorData() const
{
    if (local_to_global_index_map_->getNumberOfComponents() == 1)
    {
        // For single-variable-single-component processes reuse the existing DOF
        // table.
        const bool manage_storage = false;
        return std::make_tuple(local_to_global_index_map_.get(),
                               manage_storage);
    }

    // Otherwise construct a new DOF table.
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component;
    all_mesh_subsets_single_component.emplace_back(*mesh_subset_all_nodes_);

    const bool manage_storage = true;

    return std::make_tuple(new NumLib::LocalToGlobalIndexMap(
                               std::move(all_mesh_subsets_single_component),
                               // by location order is needed for output
                               NumLib::ComponentOrder::BY_LOCATION),
                           manage_storage);
}

void Process::initializeExtrapolator()
{
    NumLib::LocalToGlobalIndexMap* dof_table_single_component;
    bool manage_storage;

    std::tie(dof_table_single_component, manage_storage) =
        getDOFTableForExtrapolatorData();

    std::unique_ptr<NumLib::Extrapolator> extrapolator(
        new NumLib::LocalLinearLeastSquaresExtrapolator(
            *dof_table_single_component));

    // TODO: Later on the DOF table can change during the simulation!
    extrapolator_data_ = ExtrapolatorData(
        std::move(extrapolator), dof_table_single_component, manage_storage);
}

void Process::computeSparsityPattern()
{
    sparsity_pattern_ =
        NumLib::computeSparsityPattern(*local_to_global_index_map_, mesh_);
}

void Process::preTimestep(std::vector<GlobalVector*> const& x, const double t,
                          const double delta_t, const int process_id)
{
    for (auto* const solution : x)
        MathLib::LinAlg::setLocalAccessibleVector(*solution);
    preTimestepConcreteProcess(x, t, delta_t, process_id);

    boundary_conditions_[process_id].preTimestep(t, x, process_id);
}

void Process::postTimestep(std::vector<GlobalVector*> const& x, const double t,
                           const double delta_t, int const process_id)
{
    for (auto* const solution : x)
        MathLib::LinAlg::setLocalAccessibleVector(*solution);
    postTimestepConcreteProcess(x, t, delta_t, process_id);
}

void Process::postNonLinearSolver(GlobalVector const& x, const double t,
                                  double const dt, int const process_id)
{
    MathLib::LinAlg::setLocalAccessibleVector(x);
    postNonLinearSolverConcreteProcess(x, t, dt, process_id);
}

void Process::computeSecondaryVariable(double const t,
                                       double const dt,
                                       GlobalVector const& x,
                                       GlobalVector const& x_dot,
                                       int const process_id)
{
    MathLib::LinAlg::setLocalAccessibleVector(x);
    MathLib::LinAlg::setLocalAccessibleVector(x_dot);

    computeSecondaryVariableConcrete(t, dt, x, x_dot, process_id);
}

void Process::preIteration(const unsigned iter, const GlobalVector& x)
{
    MathLib::LinAlg::setLocalAccessibleVector(x);
    preIterationConcreteProcess(iter, x);
}

NumLib::IterationResult Process::postIteration(const GlobalVector& x)
{
    MathLib::LinAlg::setLocalAccessibleVector(x);
    return postIterationConcreteProcess(x);
}

}  // namespace ProcessLib
