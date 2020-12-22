/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <tuple>

#include "NumLib/ODESolver/NonlinearSolver.h"
#include "NumLib/ODESolver/ODESystem.h"
#include "NumLib/ODESolver/TimeDiscretization.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/BoundaryCondition/BoundaryConditionCollection.h"
#include "ProcessLib/Output/ExtrapolatorData.h"
#include "ProcessLib/Output/IntegrationPointWriter.h"
#include "ProcessLib/Output/SecondaryVariable.h"
#include "ProcessLib/SourceTerms/SourceTermCollection.h"

#include "AbstractJacobianAssembler.h"
#include "ProcessVariable.h"
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

    Process(std::string name_,
            MeshLib::Mesh& mesh,
            std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
            std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
                parameters,
            unsigned const integration_order,
            std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
                process_variables,
            SecondaryVariableCollection&& secondary_variables,
            const bool use_monolithic_scheme = true);

    /// Preprocessing before starting assembly for new timestep.
    void preTimestep(std::vector<GlobalVector*> const& x, const double t,
                     const double delta_t, const int process_id);

    /// Postprocessing after a complete timestep.
    void postTimestep(std::vector<GlobalVector*> const& x, const double t,
                      const double delta_t, int const process_id);

    /// Calculates secondary variables, e.g. stress and strain for deformation
    /// analysis, only after nonlinear solver being successfully conducted.
    void postNonLinearSolver(GlobalVector const& x, GlobalVector const& xdot,
                             const double t, double const dt,
                             int const process_id);

    void preIteration(const unsigned iter, GlobalVector const& x) final;

    /// compute secondary variables for the coupled equations or for output.
    void computeSecondaryVariable(double const t,
                                  double const dt,
                                  std::vector<GlobalVector*> const& x,
                                  GlobalVector const& x_dot,
                                  int const process_id);

    NumLib::IterationResult postIteration(GlobalVector const& x) final;

    void initialize();

    void setInitialConditions(
        std::vector<GlobalVector*>& process_solutions,
        std::vector<GlobalVector*>& process_solutions_prev,
        double const t,
        int const process_id);

    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int process_id) const override;

    void setCoupledSolutionsForStaggeredScheme(
        CoupledSolutionsForStaggeredScheme* const coupled_solutions)
    {
        _coupled_solutions = coupled_solutions;
    }

    void updateDeactivatedSubdomains(double const time, const int process_id);

    bool isMonolithicSchemeUsed() const { return _use_monolithic_scheme; }
    virtual void setCoupledTermForTheStaggeredSchemeToLocalAssemblers(
        int const /*process_id*/)
    {
    }

    virtual std::vector<GlobalVector> interpolateNodalValuesToIntegrationPoints(
        std::vector<GlobalVector*> const& /*nodal_values_vectors*/) const
    {
        return {};
    }

    virtual void extrapolateIntegrationPointValuesToNodes(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*integration_point_values_vectors*/,
        std::vector<GlobalVector*>& /*nodal_values_vectors*/)
    {
    }

    void preAssemble(const double t, double const dt,
                     GlobalVector const& x) final;
    void assemble(const double t, double const dt,
                  std::vector<GlobalVector*> const& x,
                  std::vector<GlobalVector*> const& xdot, int const process_id,
                  GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) final;

    void assembleWithJacobian(const double t, double const dt,
                              std::vector<GlobalVector*> const& x,
                              std::vector<GlobalVector*> const& xdot,
                              const double dxdot_dx, const double dx_dx,
                              int const process_id, GlobalMatrix& M,
                              GlobalMatrix& K, GlobalVector& b,
                              GlobalMatrix& Jac) final;

    std::vector<NumLib::IndexValueVector<GlobalIndexType>> const*
    getKnownSolutions(double const t, GlobalVector const& x,
                      int const process_id) const final
    {
        return _boundary_conditions[process_id].getKnownSolutions(t, x);
    }

    virtual NumLib::LocalToGlobalIndexMap const& getDOFTable(
        const int /*process_id*/) const
    {
        return *_local_to_global_index_map;
    }

    MeshLib::Mesh& getMesh() const { return _mesh; }
    std::vector<std::reference_wrapper<ProcessVariable>> const&
    getProcessVariables(const int process_id) const
    {
        return _process_variables[process_id];
    }

    SecondaryVariableCollection const& getSecondaryVariables() const
    {
        return _secondary_variables;
    }

    std::vector<std::unique_ptr<IntegrationPointWriter>> const*
    getIntegrationPointWriter(MeshLib::Mesh const& mesh) const
    {
        if (mesh == _mesh)
        {
            return &_integration_point_writer;
        }
        return nullptr;
    }

    // Used as a call back for CalculateSurfaceFlux process.

    virtual Eigen::Vector3d getFlux(
        std::size_t /*element_id*/,
        MathLib::Point3d const& /*p*/,
        double const /*t*/,
        std::vector<GlobalVector*> const& /*x*/) const
    {
        return Eigen::Vector3d{};
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

    /**
     * Initialize the boundary conditions for a single process or coupled
     * processes modelled by the monolithic scheme. It is called by
     * initializeBoundaryConditions().
     */
    void initializeProcessBoundaryConditionsAndSourceTerms(
        const NumLib::LocalToGlobalIndexMap& dof_table, const int process_id);

private:
    /// Process specific initialization called by initialize().
    virtual void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) = 0;

    /// Member function to initialize the boundary conditions for all coupled
    /// processes. It is called by initialize().
    virtual void initializeBoundaryConditions();

    virtual void setInitialConditionsConcreteProcess(
        std::vector<GlobalVector*>& /*x*/,
        double const /*t*/,
        int const /*process_id*/)
    {
    }

    virtual void preAssembleConcreteProcess(const double /*t*/,
                                            double const /*dt*/,
                                            GlobalVector const& /*x*/)
    {
    }

    virtual void assembleConcreteProcess(const double t, double const dt,
                                         std::vector<GlobalVector*> const& x,
                                         std::vector<GlobalVector*> const& xdot,
                                         int const process_id, GlobalMatrix& M,
                                         GlobalMatrix& K, GlobalVector& b) = 0;

    virtual void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
        const double dx_dx, int const process_id, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) = 0;

    virtual void preTimestepConcreteProcess(
        std::vector<GlobalVector*> const& /*x*/,
        const double /*t*/,
        const double /*dt*/,
        const int /*process_id*/)
    {
    }

    virtual void postTimestepConcreteProcess(
        std::vector<GlobalVector*> const& /*x*/,
        const double /*t*/,
        const double /*dt*/,
        int const /*process_id*/)
    {
    }

    virtual void postNonLinearSolverConcreteProcess(
        GlobalVector const& /*x*/, GlobalVector const& /*xdot*/,
        const double /*t*/, double const /*dt*/, int const /*process_id*/)
    {
    }

    virtual void preIterationConcreteProcess(const unsigned /*iter*/,
                                             GlobalVector const& /*x*/)
    {
    }

    virtual void computeSecondaryVariableConcrete(
        double const /*t*/,
        double const /*dt*/,
        std::vector<GlobalVector*> const& /*x*/,
        GlobalVector const& /*x_dot*/,
        int const /*process_id*/)
    {
    }

    virtual NumLib::IterationResult postIterationConcreteProcess(
        GlobalVector const& /*x*/)
    {
        return NumLib::IterationResult::SUCCESS;
    }

protected:
    /** This function is for general cases, in which all equations of the
     coupled processes have the same number of unknowns. For the general cases
     with the staggered scheme, all equations of the coupled processes share one
     DOF table hold by  @c _local_to_global_index_map. Other cases can be
     considered by overloading this member function in the derived class.
     */
    virtual void constructDofTable();

    /**
     * Construct the DOF table for the monolithic scheme,
     * which is stored in the
     * member of this class, @c _local_to_global_index_map.
     */
    void constructMonolithicProcessDofTable();

    /**
     * Construct the DOF table for a specified process in the staggered scheme,
     * which is stored in the
     * member of this class, @c _local_to_global_index_map.
     */
    void constructDofTableOfSpecifiedProsessStaggeredScheme(
        const int specified_prosess_id);

    /**
     * Get the address of a LocalToGlobalIndexMap, and the status of its memory.
     * If the LocalToGlobalIndexMap is created as new in this function, the
     * function also returns a true boolean value to let Extrapolator manage
     * the memory by the address returned by this function.
     *
     * @return Address of a LocalToGlobalIndexMap and its memory status.
     */
    virtual std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
    getDOFTableForExtrapolatorData() const;

private:
    void initializeExtrapolator();

    /// Computes and stores global matrix' sparsity pattern from given
    /// DOF-table.
    void computeSparsityPattern();

public:
    std::string const name;

protected:
    MeshLib::Mesh& _mesh;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_all_nodes;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _local_to_global_index_map;

    SecondaryVariableCollection _secondary_variables;

    VectorMatrixAssembler _global_assembler;

    const bool _use_monolithic_scheme;

    /// Pointer to CoupledSolutionsForStaggeredScheme, which contains the
    /// references to the solutions of the coupled processes.
    CoupledSolutionsForStaggeredScheme* _coupled_solutions;

    /// Order of the integration method for element-wise integration.
    /// The Gauss-Legendre integration method and available orders is
    /// implemented in MathLib::GaussLegendre.
    unsigned const _integration_order;

    /// An optional vector containing descriptions for integration point data
    /// output and setting of the integration point initial conditions.
    /// The integration point writer are implemented in specific processes.
    std::vector<std::unique_ptr<IntegrationPointWriter>>
        _integration_point_writer;

    GlobalSparsityPattern _sparsity_pattern;

protected:
    /// Variables used by this process.  For the monolithic scheme or a
    /// single process, the size of the outer vector is one. For the
    /// staggered scheme, the size of the outer vector is the number of the
    /// coupled processes.
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        _process_variables;

    /// Vector for boundary conditions. For the monolithic scheme or a
    /// single process, the size of the vector is one. For the staggered
    /// scheme, the size of vector is the number of the coupled processes.
    std::vector<BoundaryConditionCollection> _boundary_conditions;

private:
    /// Vector for nodal source term collections. For the monolithic scheme
    /// or a single process, the size of the vector is one. For the staggered
    /// scheme, the size of vector is the number of the coupled processes.
    std::vector<SourceTermCollection> _source_term_collections;

    ExtrapolatorData _extrapolator_data;
};

}  // namespace ProcessLib
