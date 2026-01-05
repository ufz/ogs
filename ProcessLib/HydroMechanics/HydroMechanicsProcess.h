// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "HydroMechanicsProcessData.h"
#include "LocalAssemblerInterface.h"
#include "ProcessLib/AssemblyMixin.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace HydroMechanics
{
/// Linear kinematics poro-mechanical/biphasic (fluid-solid mixture) model.
///
/// The mixture momentum balance and the mixture mass balance are solved under
/// fully saturated conditions.
template <int DisplacementDim>
class HydroMechanicsProcess final
    : public Process,
      private AssemblyMixin<HydroMechanicsProcess<DisplacementDim>>
{
    friend class AssemblyMixin<HydroMechanicsProcess<DisplacementDim>>;

public:
    HydroMechanicsProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        HydroMechanicsProcessData<DisplacementDim>&& process_data,
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
     *                               process_id = 0. For the staggered scheme,
     *                               process_id = 0 represents the
     *                               hydraulic (H) process, while process_id = 1
     *                               represents the mechanical (M) process.
     * @return Matrix specifications including size and sparse pattern.
     */
    MathLib::MatrixSpecifications getMatrixSpecifications(
        const int process_id) const override;

private:
    using LocalAssemblerIF = LocalAssemblerInterface<DisplacementDim>;

    void constructDofTable() override;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void initializeBoundaryConditions(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media) override;

    void assembleConcreteProcess(const double t, double const /*dt*/,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& x_prev,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const /*dt*/,
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalVector& b, GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                    double const t, double const dt,
                                    const int process_id) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     std::vector<GlobalVector*> const& x_prev,
                                     const double t, const double dt,
                                     int const process_id) override;

    void postNonLinearSolverConcreteProcess(
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, const double t,
        double const dt, int const process_id) override;

    void setInitialConditionsConcreteProcess(std::vector<GlobalVector*>& x,
                                             double const t,
                                             int const process_id) override;

    NumLib::LocalToGlobalIndexMap const& getDOFTable(
        const int process_id) const override;

    std::vector<std::vector<std::string>> initializeAssemblyOnSubmeshes(
        std::vector<std::reference_wrapper<MeshLib::Mesh>> const& meshes)
        override;

    bool isMonolithicSchemeUsed() const override
    {
        return process_data_.isMonolithicSchemeUsed();
    }

private:
    std::vector<MeshLib::Node*> base_nodes_;
    std::unique_ptr<MeshLib::MeshSubset const> mesh_subset_base_nodes_;
    HydroMechanicsProcessData<DisplacementDim> process_data_;

    std::vector<std::unique_ptr<LocalAssemblerIF>> local_assemblers_;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        local_to_global_index_map_single_component_;

    /// Local to global index mapping for base nodes, which is used for linear
    /// interpolation for pressure in the staggered scheme.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        local_to_global_index_map_with_base_nodes_;

    /// Sparsity pattern for the flow equation, and it is initialized only if
    /// the staggered scheme is used.
    GlobalSparsityPattern sparsity_pattern_with_linear_element_;

    void computeSecondaryVariableConcrete(double const t, double const dt,
                                          std::vector<GlobalVector*> const& x,
                                          GlobalVector const& x_prev,
                                          const int process_id) override;
    /**
     * @copydoc ProcessLib::Process::getDOFTableForExtrapolatorData()
     */
    std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
    getDOFTableForExtrapolatorData() const override;

    /// Check whether the process represented by \c process_id is/has
    /// mechanical process. In the present implementation, the mechanical
    /// process has process_id == 1 in the staggered scheme.
    bool hasMechanicalProcess(int const process_id) const
    {
        return process_id == process_data_.mechanics_related_process_id;
    }
};

extern template class HydroMechanicsProcess<2>;
extern template class HydroMechanicsProcess<3>;

}  // namespace HydroMechanics
}  // namespace ProcessLib
