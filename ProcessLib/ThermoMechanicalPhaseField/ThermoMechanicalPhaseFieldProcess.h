/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Process.h"
#include "LocalAssemblerInterface.h"

#include "ThermoMechanicalPhaseFieldFEM.h"
#include "ThermoMechanicalPhaseFieldProcessData.h"

namespace ProcessLib
{
namespace ThermoMechanicalPhaseField
{
struct ThermoMechanicalPhaseFieldLocalAssemblerInterface;

template <int DisplacementDim>
class ThermoMechanicalPhaseFieldProcess final : public Process
{
    using Base = Process;

public:
    ThermoMechanicalPhaseFieldProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        ThermoMechanicalPhaseFieldProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller,
        bool const use_monolithic_scheme,
        int const mechanics_related_process_id,
        int const phase_field_process_id,
        int const heat_conduction_process_id);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;
    //! @}

    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int process_id) const override;

    NumLib::LocalToGlobalIndexMap const& getDOFTable(
        const int process_id) const override;

private:
    void constructDofTable() override;

    void initializeBoundaryConditions() override;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                                    double const dt,
                                    const int process_id) override;

    void postTimestepConcreteProcess(GlobalVector const& x,
                                     int const process_id) override;

    void postNonLinearSolverConcreteProcess(GlobalVector const& x,
                                            const double t,
                                            int const process_id) override;

    // To be replaced.
    NumLib::LocalToGlobalIndexMap& getDOFTableByProcessID(
        const int process_id) const;

private:
    ThermoMechanicalPhaseFieldProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<ThermoMechanicalPhaseFieldLocalAssemblerInterface>>
        _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;

    /// Sparsity pattern for the phase field equation, and it is initialized
    //  only if the staggered scheme is used.
    GlobalSparsityPattern _sparsity_pattern_with_single_component;

    /// ID of the processes that contains mechanical process.
    int const _mechanics_related_process_id;

    /// ID of phase field process.
    int const _phase_field_process_id;

    /// ID of heat conduction process.
    int const _heat_conduction_process_id;
};

extern template class ThermoMechanicalPhaseFieldProcess<2>;
extern template class ThermoMechanicalPhaseFieldProcess<3>;

}  // namespace ThermoMechanicalPhaseField
}  // namespace ProcessLib
