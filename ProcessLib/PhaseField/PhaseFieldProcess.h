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
#include "PhaseFieldProcessData.h"

namespace ProcessLib
{
namespace PhaseField
{
/**
 * \brief A class to simulate mechanical fracturing process
 * using phase-field approach in solids described by
 *
 * \f[
 *     \mathrm{div} \left[ \left(d^2 + k \right) \boldsymbol{\sigma}_0^+
 *      + \boldsymbol{\sigma}_0^- \right] +
 *        \varrho \boldsymbol{b} = \boldsymbol{0}
 * \f]
 * \f[
 *     2d H(\boldsymbol{\epsilon}_\mathrm{el})
 *      - \frac{1 - d}{2 \varepsilon} g_\mathrm{c}
 *      - 2 \varepsilon g_\mathrm{c} \mathrm{div}(\mathrm{grad} d) = 0
 * \f]
 * where
 *    \f{eqnarray*}{
 *       &d:&                    \mbox{order parameter,}\\
 *       &\varrho:&              \mbox{density,}\\
 *       &H:&                    \mbox{history field,}\\
 *       &g_\mathrm{c}:&         \mbox{fracture energy,}\\
 *       &\varepsilon:&          \mbox{length scale}\\
 *    \f}
 *
 * Detailed model description can refer
 * <a href="Miao_Biot2017.pdff" target="_blank"><b>Phase field method</b></a>
 */
template <int DisplacementDim>
class PhaseFieldProcess final : public Process
{
    using Base = Process;

public:
    PhaseFieldProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        PhaseFieldProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller,
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

private:
    PhaseFieldProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;

    MeshLib::PropertyVector<double>* _nodal_forces = nullptr;

    /// Sparsity pattern for the phase field equation, and it is initialized
    ///  only if the staggered scheme is used.
    GlobalSparsityPattern _sparsity_pattern_with_single_component;

    /// Check whether the process represented by \c process_id is/has
    /// mechanical process. In the present implementation, the mechanical
    /// process has process_id == 0 in the staggered scheme.
    bool isPhaseFieldProcess(int const process_id) const
    {
        return !_use_monolithic_scheme && process_id == 1;
    }
};

extern template class PhaseFieldProcess<2>;
extern template class PhaseFieldProcess<3>;

}  // namespace PhaseField
}  // namespace ProcessLib
