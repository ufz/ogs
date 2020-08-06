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
#include "HydroMechanicsProcessData.h"
#include "LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace HydroMechanics
{
/// Linear kinematics poro-mechanical/biphasic (fluid-solid mixture) model.
///
/// The mixture momentum balance and the mixture mass balance are solved under
/// fully saturated conditions.
template <int DisplacementDim>
class HydroMechanicsProcess final : public Process
{
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

    void initializeBoundaryConditions() override;

    void assembleConcreteProcess(const double t, double const /*dt*/,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& xdot,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const /*dt*/,
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
        const double dx_dx, int const process_id, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                    double const t, double const dt,
                                    const int process_id) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     const double t, const double dt,
                                     int const process_id) override;

    void postNonLinearSolverConcreteProcess(GlobalVector const& x,
                                            const double t, double const dt,
                                            int const process_id) override;

    NumLib::LocalToGlobalIndexMap const& getDOFTable(
        const int process_id) const override;

private:
    std::vector<MeshLib::Node*> _base_nodes;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_base_nodes;
    HydroMechanicsProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerIF>> _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;

    /// Local to global index mapping for base nodes, which is used for linear
    /// interpolation for pressure in the staggered scheme.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_with_base_nodes;

    /// Sparsity pattern for the flow equation, and it is initialized only if
    /// the staggered scheme is used.
    GlobalSparsityPattern _sparsity_pattern_with_linear_element;

    /// Solutions of the previous time step
    std::array<std::unique_ptr<GlobalVector>, 2> _xs_previous_timestep;

    void computeSecondaryVariableConcrete(double const t, double const dt,
                                          GlobalVector const& x,
                                          GlobalVector const& x_dot,
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
        return _use_monolithic_scheme || process_id == 1;
    }

    MeshLib::PropertyVector<double>* _nodal_forces = nullptr;
    MeshLib::PropertyVector<double>* _hydraulic_flow = nullptr;
};

extern template class HydroMechanicsProcess<2>;
extern template class HydroMechanicsProcess<3>;

}  // namespace HydroMechanics
}  // namespace ProcessLib
