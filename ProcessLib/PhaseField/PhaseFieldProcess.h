/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "LocalAssemblerInterface.h"
#include "PhaseFieldProcessData.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace PhaseField
{
template <int DisplacementDim>
class PhaseFieldProcess final : public Process
{
public:
    PhaseFieldProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        PhaseFieldProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        bool const use_monolithic_scheme);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;
    //! @}

    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int process_id) const override;

    NumLib::LocalToGlobalIndexMap const& getDOFTable(
        const int process_id) const override;

private:
    using LocalAssemblerInterface = PhaseFieldLocalAssemblerInterface;

    void constructDofTable() override;

    void initializeBoundaryConditions(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media) override;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& x_prev,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalVector& b, GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                    double const t, double const dt,
                                    const int process_id) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     std::vector<GlobalVector*> const& x_prev,
                                     const double t, const double delta_t,
                                     int const process_id) override;

    void postNonLinearSolverConcreteProcess(
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, const double t,
        double const dt, int const process_id) override;

    void updateConstraints(GlobalVector& lower, GlobalVector& upper,
                           int const process_id) override;

private:
    PhaseFieldProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;

    MeshLib::PropertyVector<double>* _nodal_forces = nullptr;

    /// Sparsity pattern for the phase field equation, and it is initialized
    ///  only if the staggered scheme is used.
    GlobalSparsityPattern _sparsity_pattern_with_single_component;

    std::unique_ptr<GlobalVector> _x_previous_timestep;

    /// Check whether the process represented by \c process_id is/has
    /// mechanical process. In the present implementation, the mechanical
    /// process has process_id == 0 in the staggered scheme.
    bool isPhaseFieldProcess(int const process_id) const;
};

extern template class PhaseFieldProcess<2>;
extern template class PhaseFieldProcess<3>;

}  // namespace PhaseField
}  // namespace ProcessLib
