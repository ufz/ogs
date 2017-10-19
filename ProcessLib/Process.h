/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/ODESolver/NonlinearSolver.h"
#include "NumLib/ODESolver/ODESystem.h"
#include "NumLib/ODESolver/TimeDiscretization.h"
#include "NumLib/NamedFunctionCaller.h"
#include "ProcessLib/BoundaryCondition/BoundaryConditionCollection.h"
#include "ProcessLib/SourceTerms/NodalSourceTerm.h"
#include "ProcessLib/Parameter/Parameter.h"

#include "ExtrapolatorData.h"
#include "ProcessVariable.h"
#include "SecondaryVariable.h"
#include "CachedSecondaryVariable.h"
#include "AbstractJacobianAssembler.h"
#include "VectorMatrixAssembler.h"

namespace MeshLib
{
class Mesh;
}

namespace ProcessLib
{
struct CoupledSolutionsForStaggeredScheme;

class Process
    : public NumLib::ODESystem<  // TODO: later on use a simpler ODE system
          NumLib::ODESystemTag::FirstOrderImplicitQuasilinear,
          NumLib::NonlinearSolverTag::Newton>
{
public:
    using NonlinearSolver = NumLib::NonlinearSolverBase;
    using TimeDiscretization = NumLib::TimeDiscretization;

    Process(MeshLib::Mesh& mesh,
            std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
            std::vector<std::unique_ptr<ParameterBase>> const& parameters,
            unsigned const integration_order,
            std::vector<std::reference_wrapper<ProcessVariable>>&&
                process_variables,
            SecondaryVariableCollection&& secondary_variables,
            NumLib::NamedFunctionCaller&& named_function_caller);

    /// Preprocessing before starting assembly for new timestep.
    void preTimestep(GlobalVector const& x, const double t,
                     const double delta_t, const int process_id);

    /// Postprocessing after a complete timestep.
    void postTimestep(GlobalVector const& x);

    void preIteration(const unsigned iter, GlobalVector const& x) final;

    /// compute secondary variables for the coupled equations or for output.
    void computeSecondaryVariable(const double t, GlobalVector const& x);

    NumLib::IterationResult postIteration(GlobalVector const& x) final;

    void initialize();

    void setInitialConditions(const unsigned pcs_id, const double t,
                              GlobalVector& x);

    MathLib::MatrixSpecifications getMatrixSpecifications() const final;

    void setCoupledSolutionsForStaggeredScheme(
        CoupledSolutionsForStaggeredScheme* const coupled_solutions)
    {
        _coupled_solutions = coupled_solutions;

    }
    void setDecouplingSchemeType(const bool is_monolithic_scheme)
    {
        _is_monolithic_scheme = is_monolithic_scheme;
    }

    bool useMonolithicScheme() const { return _is_monolithic_scheme; }

    virtual void setCoupledSolutionsForStaggeredSchemeToLocalAssemblers() {}
    void preAssemble(const double t, GlobalVector const& x) override final;
    void assemble(const double t, GlobalVector const& x, GlobalMatrix& M,
                  GlobalMatrix& K, GlobalVector& b) final;

    void assembleWithJacobian(const double t, GlobalVector const& x,
                              GlobalVector const& xdot, const double dxdot_dx,
                              const double dx_dx, GlobalMatrix& M,
                              GlobalMatrix& K, GlobalVector& b,
                              GlobalMatrix& Jac) final;

    std::vector<NumLib::IndexValueVector<GlobalIndexType>> const*
    getKnownSolutions(double const t) const final
    {
        return _boundary_conditions.getKnownSolutions(t);
    }

    NumLib::LocalToGlobalIndexMap const& getDOFTable() const
    {
        return *_local_to_global_index_map;
    }

    MeshLib::Mesh& getMesh() const { return _mesh; }
    std::vector<std::reference_wrapper<ProcessVariable>> const&
    getProcessVariables() const
    {
        return _process_variables;
    }

    SecondaryVariableCollection const& getSecondaryVariables() const
    {
        return _secondary_variables;
    }

    // Get the solution of the previous time step.
    virtual GlobalVector* getPreviousTimeStepSolution(
        const int /*process_id*/) const
    {
        return nullptr;
    }

    // Used as a call back for CalculateSurfaceFlux process.
    virtual std::vector<double> getFlux(std::size_t /*element_id*/,
                                        MathLib::Point3d const& /*p*/,
                                        GlobalVector const& /*x*/) const
    {
        return std::vector<double>{};
    }

protected:
    NumLib::Extrapolator& getExtrapolator() const
    {
        return _extrapolator_data.getExtrapolator();
    }

    NumLib::LocalToGlobalIndexMap const& getSingleComponentDOFTable() const
    {
        return _extrapolator_data.getDOFTable();
    }

private:
    /// Process specific initialization called by initialize().
    virtual void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) = 0;

    virtual void preAssembleConcreteProcess(const double /*t*/,
                                            GlobalVector const& /*x*/)
    {
    }

    virtual void assembleConcreteProcess(const double t, GlobalVector const& x,
                                         GlobalMatrix& M, GlobalMatrix& K,
                                         GlobalVector& b) = 0;

    virtual void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) = 0;

    virtual void preTimestepConcreteProcess(GlobalVector const& /*x*/,
                                            const double /*t*/,
                                            const double /*delta_t*/,
                                            const int /*process_id*/)
    {
    }

    virtual void postTimestepConcreteProcess(GlobalVector const& /*x*/) {}
    virtual void preIterationConcreteProcess(const unsigned /*iter*/,
                                             GlobalVector const& /*x*/)
    {
    }

    virtual void computeSecondaryVariableConcrete(const double /*t*/,
                                                  GlobalVector const& /*x*/)
    {
    }

    virtual NumLib::IterationResult postIterationConcreteProcess(
        GlobalVector const& /*x*/)
    {
        return NumLib::IterationResult::SUCCESS;
    }

protected:
    virtual void constructDofTable();

private:
    void initializeExtrapolator();

    /// Finishes the \c _named_function_caller and adds a secondary variable for
    /// each of the named functions.
    void finishNamedFunctionsInitialization();

    /// Computes and stores global matrix' sparsity pattern from given
    /// DOF-table.
    void computeSparsityPattern();

    void setVariableInitialCondition(const int variable_id, double const t,
                                     GlobalVector& x);

protected:
    MeshLib::Mesh& _mesh;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_all_nodes;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _local_to_global_index_map;

    SecondaryVariableCollection _secondary_variables;

    NumLib::NamedFunctionCaller _named_function_caller;
    std::vector<std::unique_ptr<CachedSecondaryVariable>>
        _cached_secondary_variables;
    SecondaryVariableContext _secondary_variable_context;

    VectorMatrixAssembler _global_assembler;

    mutable bool _is_monolithic_scheme;

    /// Pointer to CoupledSolutionsForStaggeredScheme, which contains the
    /// references to the solutions of the coupled processes.
    CoupledSolutionsForStaggeredScheme* _coupled_solutions;

    /// Set the solutions of the previous time step to the coupled term.
    /// It only performs for the staggered scheme.
    void setCoupledSolutionsOfPreviousTimeStep();

    /// Order of the integration method for element-wise integration.
    /// The Gauss-Legendre integration method and available orders is
    /// implemented in MathLib::GaussLegendre.
    unsigned const _integration_order;

private:
    GlobalSparsityPattern _sparsity_pattern;

    /// Variables used by this process.
    std::vector<std::reference_wrapper<ProcessVariable>> _process_variables;

    BoundaryConditionCollection _boundary_conditions;
    std::vector<std::unique_ptr<NodalSourceTerm>> _source_terms;

    ExtrapolatorData _extrapolator_data;
};

}  // namespace ProcessLib
