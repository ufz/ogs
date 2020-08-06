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

#include "ProcessLib/Process.h"

#include "LocalAssemblerInterface.h"
#include "ThermoMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoMechanics
{
template <int DisplacementDim>
class ThermoMechanicsProcess final : public Process
{
public:
    ThermoMechanicsProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        ThermoMechanicsProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        bool const use_monolithic_scheme);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;
    //! @}

    /**
     * Get the size and the sparse pattern of the global matrix in order to
     * create the global matrices and vectors for the system equations of this
     * process.
     *
     * @param process_id Process ID. If the monolithic scheme is applied,
     *                               process_id = 0.
     * @return Matrix specifications including size and sparse pattern.
     */
    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int process_id) const override;

private:
    using LocalAssemblerInterface =
        ThermoMechanicsLocalAssemblerInterface<DisplacementDim>;

    void constructDofTable() override;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void initializeBoundaryConditions() override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& xdot,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
        const double dx_dx, int const process_id, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                    double const t, double const dt,
                                    const int process_id) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     const double t, const double dt,
                                     int const process_id) override;

    NumLib::LocalToGlobalIndexMap const& getDOFTable(
        const int process_id) const override;

private:
    ThermoMechanicsProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;

    /// Sparsity pattern for the heat conduction equation, and it is initialized
    /// only if the staggered scheme is used.
    GlobalSparsityPattern _sparsity_pattern_with_single_component;

    MeshLib::PropertyVector<double>* _nodal_forces = nullptr;
    MeshLib::PropertyVector<double>* _heat_flux = nullptr;

    /// Temperature of the previous time step for staggered scheme.
    std::unique_ptr<GlobalVector> _previous_T;

    /// Set the increment solutions of the present time step to the coupled
    /// term.
    void setCoupledSolutionsOfPreviousTimeStep();
};

extern template class ThermoMechanicsProcess<2>;
extern template class ThermoMechanicsProcess<3>;

}  // namespace ThermoMechanics
}  // namespace ProcessLib
