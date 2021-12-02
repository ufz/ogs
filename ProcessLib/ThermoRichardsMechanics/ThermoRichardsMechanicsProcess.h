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

#include "LocalAssemblerInterface.h"
#include "ProcessLib/Process.h"
#include "ThermoRichardsMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoRichardsMechanics
{
/**
 * \brief Global assembler for the monolithic scheme of the non-isothermal
 * Richards flow coupled with mechanics.
 *
 * <b>Governing equations without vapor diffusion</b>
 *
 * The energy balance equation is given by
 * \f[
 *  (\rho c_p)^{eff}\dot T -
 *  \nabla (\mathbf{k}_T^{eff} \nabla T)+\rho^l c_p^l \nabla T \cdot
 * \mathbf{v}^l
 * = Q_T
 * \f]
 *  with\f$T\f$ the temperature, \f$(\rho c_p)^{eff}\f$ the  effective
 * volumetric heat
 * capacity, \f$\mathbf{k}_T^{eff} \f$
 *  the effective thermal conductivity, \f$\rho^l\f$ the density of liquid,
 * \f$c_p^l\f$ the specific heat  capacity of liquid, \f$\mathbf{v}^l\f$ the
 * liquid velocity, and \f$Q_T\f$ the point heat source. The  effective
 * volumetric heat can be considered as a composite of the contributions of
 * solid phase and the liquid phase as
 * \f[
 * (\rho c_p)^{eff} = (1-\phi) \rho^s c_p^s + S^l \phi \rho^l c_p^l
 * \f]
 * with \f$\phi\f$ the porosity, \f$S^l\f$  the liquid saturation, \f$\rho^s \f$
 * the solid density, and \f$c_p^s\f$ the specific heat capacity of solid.
 * Similarly, the effective thermal conductivity is given by
 * \f[
 * \mathbf{k}_T^{eff} = (1-\phi) \mathbf{k}_T^s + S^l \phi k_T^l \mathbf I
 * \f]
 * where \f$\mathbf{k}_T^s\f$ is the thermal conductivity tensor of solid, \f$
 *  k_T^l\f$ is the thermal conductivity of liquid, and \f$\mathbf I\f$ is the
 * identity tensor.
 *
 * The mass balance equation is given by
 * \f{eqnarray*}{
 * \left(S^l\beta - \phi\frac{\partial S}{\partial p_c}\right) \rho^l\dot p
 * - S \left( \frac{\partial \rho^l}{\partial T}
 * +\rho^l(\alpha_B -S)
 * \alpha_T^s
 * \right)\dot T\\
 *  +\nabla (\rho^l \mathbf{v}^l) + S \alpha_B \rho^l \nabla \cdot \dot {\mathbf
 * u}= Q_H
 * \f}
 * where \f$p\f$ is the pore pressure,  \f$p_c\f$ is the
 * capillary pressure, which is \f$-p\f$ under the single phase assumption,
 *  \f$\beta\f$ is a composite coefficient by the liquid compressibility and
 * solid compressibility, \f$\alpha_B\f$ is the Biot's constant,
 * \f$\alpha_T^s\f$ is the linear thermal  expansivity of solid, \f$Q_H\f$
 * is the point source or sink term,  \f$ \mathbf u\f$ is the displacement, and
 * \f$H(S-1)\f$ is the Heaviside function.
 * The liquid velocity \f$\mathbf{v}^l\f$ is
 * described by the Darcy's law as
 * \f[
 * \mathbf{v}^l=-\frac{{\mathbf k} k_{ref}}{\mu} (\nabla p - \rho^l \mathbf g)
 * \f]
 * with \f${\mathbf k}\f$ the intrinsic permeability, \f$k_{ref}\f$ the relative
 * permeability, \f$\mathbf g\f$ the gravitational force.
 *
 * The momentum balance equation takes the form of
 * \f[
 * \nabla (\mathbf{\sigma}-b(S)\alpha_B p^l \mathbf I) +\mathbf f=0
 * \f]
 * with \f$\mathbf{\sigma}\f$  the effective stress tensor, \f$b(S)\f$ the
 * Bishop model, \f$\mathbf f\f$ the body force, and \f$\mathbf I\f$ the
 * identity. The primary unknowns of the momentum balance equation are the
 * displacement \f$\mathbf u\f$, which is associated with the stress by the
 * the generalized Hook's law as
 * \f[
 * {\dot {\mathbf {\sigma}}} = C {\dot {\mathbf \epsilon}}^e
 *  = C ( {\dot {\mathbf \epsilon}} - {\dot {\mathbf \epsilon}}^T
 * -{\dot {\mathbf \epsilon}}^p - {\dot {\mathbf \epsilon}}^{sw}-\cdots )
 * \f]
 *  with \f$C\f$ the forth order elastic tensor,
 *  \f${\dot {\mathbf \epsilon}}\f$ the total strain rate,
 *  \f${\dot {\mathbf \epsilon}}^e\f$ the elastic strain rate,
 *  \f${\dot {\mathbf \epsilon}}^T\f$ the thermal strain rate,
 *  \f${\dot {\mathbf \epsilon}}^p\f$ the plastic strain rate,
 *  \f${\dot {\mathbf \epsilon}}^{sw}\f$ the swelling strain rate.
 *
 *  The strain tensor is given by displacement vector as
 *  \f[
 *   \mathbf \epsilon =
 *   \frac{1}{2} \left((\nabla \mathbf u)^{\text T}+\nabla \mathbf u\right)
 * \f]
 * where the superscript \f${\text T}\f$ means transpose,
 */
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
    using LocalAssemblerIF = LocalAssemblerInterface<DisplacementDim>;

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
        std::vector<GlobalVector*> const& xdot, int const process_id,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac) override;

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
    MeshLib::PropertyVector<double>* heat_flux_ = nullptr;
};

extern template class ThermoRichardsMechanicsProcess<2>;
extern template class ThermoRichardsMechanicsProcess<3>;

}  // namespace ThermoRichardsMechanics
}  // namespace ProcessLib
