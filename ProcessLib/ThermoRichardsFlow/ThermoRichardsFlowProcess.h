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
#include "ThermoRichardsFlowProcessData.h"

namespace ProcessLib
{
namespace ThermoRichardsFlow
{
/**
 * \brief Global assembler for the monolithic scheme of the non-isothermal
 * Richards flow.
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
 * is the point source or sink term,  \f$H(S-1)\f$ is the Heaviside function, and
 * \f$ \mathbf u\f$ is the displacement. While this process does not contain a fully
 * mechanical coupling, simplfied expressions can be given to approximate the latter
 * term under certain stress conditions.
 * The liquid velocity \f$\mathbf{v}^l\f$ is
 * described by the Darcy's law as
 * \f[
 * \mathbf{v}^l=-\frac{{\mathbf k} k_{ref}}{\mu} (\nabla p - \rho^l \mathbf g)
 * \f]
 * with \f${\mathbf k}\f$ the intrinsic permeability, \f$k_{ref}\f$ the relative
 * permeability, \f$\mathbf g\f$ the gravitational force.
 */
class ThermoRichardsFlowProcess final : public Process
{
public:
    ThermoRichardsFlowProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        ThermoRichardsFlowProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        bool const use_monolithic_scheme);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

private:
    using LocalAssemblerIF = LocalAssemblerInterface;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void setInitialConditionsConcreteProcess(std::vector<GlobalVector*>& x,
                                             double const t,
                                             int const /*process_id*/) override;

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

private:
    std::vector<MeshLib::Node*> _base_nodes;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_base_nodes;
    ThermoRichardsFlowProcessData _process_data;

    std::vector<std::unique_ptr<LocalAssemblerIF>> _local_assemblers;

    /// Sparsity pattern for the flow equation, and it is initialized only if
    /// the staggered scheme is used.
    GlobalSparsityPattern _sparsity_pattern_with_linear_element;

    void computeSecondaryVariableConcrete(double const t, double const dt,
                                          std::vector<GlobalVector*> const& x,
                                          GlobalVector const& x_dot,
                                          int const process_id) override;
    std::vector<NumLib::LocalToGlobalIndexMap const*> getDOFTables(
        const int number_of_processes) const;

    MeshLib::PropertyVector<double>* _heat_flux = nullptr;
    MeshLib::PropertyVector<double>* _hydraulic_flow = nullptr;
};

}  // namespace ThermoRichardsFlow
}  // namespace ProcessLib
