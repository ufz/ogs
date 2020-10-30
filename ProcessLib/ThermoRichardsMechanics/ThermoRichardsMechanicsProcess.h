/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Process.h"
#include "ProcessLib/RichardsMechanics/LocalAssemblerInterface.h"
#include "ThermoRichardsMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoRichardsMechanics
{
/// Global assembler for the monolithic scheme of the non-isothermal Richards
/// flow coupled with mechanics.
///
template <int DisplacementDim>
class ThermoRichardsMechanicsProcess final : public Process
{
public:
    ThermoRichardsMechanicsProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        ThermoRichardsMechanicsProcessData<DisplacementDim>&& process_data,
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
    using LocalAssemblerIF =
        RichardsMechanics::LocalAssemblerInterface<DisplacementDim>;

    void constructDofTable() override;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void initializeBoundaryConditions() override;

    void setInitialConditionsConcreteProcess(std::vector<GlobalVector*>& x,
                                             double const t,
                                             int const process_id) override;

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

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     double const t, double const dt,
                                     const int process_id) override;

    NumLib::LocalToGlobalIndexMap const& getDOFTable(
        const int process_id) const override;

private:
    std::vector<MeshLib::Node*> base_nodes_;
    std::unique_ptr<MeshLib::MeshSubset const> mesh_subset_base_nodes_;
    ThermoRichardsMechanicsProcessData<DisplacementDim> process_data_;

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
                                          GlobalVector const& x_dot,
                                          int const process_id) override;
    /**
     * @copydoc ProcessLib::Process::getDOFTableForExtrapolatorData()
     */
    std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
    getDOFTableForExtrapolatorData() const override;

    std::vector<NumLib::LocalToGlobalIndexMap const*> getDOFTables(
        const int number_of_processes) const;

    /// Check whether the process represented by \c process_id is/has
    /// mechanical process. In the present implementation, the mechanical
    /// process has process_id == 1 in the staggered scheme.
    bool hasMechanicalProcess(int const /*process_id*/) const { return true; }

    MeshLib::PropertyVector<double>* nodal_forces_ = nullptr;
    MeshLib::PropertyVector<double>* hydraulic_flow_ = nullptr;
};

extern template class ThermoRichardsMechanicsProcess<2>;
extern template class ThermoRichardsMechanicsProcess<3>;

}  // namespace ThermoRichardsMechanics
}  // namespace ProcessLib
