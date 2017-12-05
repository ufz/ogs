/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "RichardsComponentTransportFEM.h"
#include "RichardsComponentTransportProcessData.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace RichardsComponentTransport
{
/**
 * # RichardsComponentTransport process
 *
 * ## Governing Equations
 *
 * ### Richards Flow
 *
 * The flow process is described by
 * \f[
 * \phi \frac{\partial \rho_w}{\partial p} \frac{\partial p}{\partial t} S
 *     - \phi \rho_w \frac{\partial S}{\partial p_c}
 *       \frac{\partial p_c}{\partial t}
 *     - \nabla \cdot \left[\rho_w \frac{k_{\mathrm{rel}} \kappa}{\mu}
 *       \nabla \left( p + \rho_w g z \right)\right]
 *     - Q_p = 0,
 * \f]
 * where
 * - \f$\phi\f$ is the porosity,
 * - \f$S\f$ is the saturation,
 * - \f$p\f$ is the pressure,
 * - \f$k_{\mathrm{rel}}\f$ is the relativ permeability (depending on \f$S\f$),
 * - \f$\kappa\f$ is the specific permeability,
 * - \f$\mu\f$ is viscosity of the fluid,
 * - \f$\rho_w\f$ is the mass density of the fluid, and
 * - \f$g\f$ is the gravitational acceleration.
 *
 * Here it is assumed, that
 * - the porosity is constant and
 * - the pressure of the gas phase is zero.
 *
 * The capillary pressure is given by
 * \f[
 * p_c = \frac{\rho_w g}{\alpha}
 *  \left[S_{\mathrm{eff}}^{-\frac{1}{m}} - 1\right]^{\frac{1}{n}}
 * \f]
 * and the effective saturation by
 * \f[
 * S_{\mathrm{eff}} = \frac{S-S_r}{S_{\max} - S_r}
 * \f]
 *
 * ### Mass Transport
 * The mass transport process is described by
 * \f[
 * \phi R \frac{\partial C}{\partial t}
    + \nabla \cdot \left(\vec{q} C - D \nabla C \right)
        + \phi R \vartheta C - Q_C = 0
 * \f]
 * where
 * - \f$R\f$ is the retardation factor,
 * - \f$C\f$ is the concentration,
 * - \f$\vec{q} = \frac{k_{\mathrm{rel}} \kappa}{\mu(C)}
 *   \nabla \left( p + \rho_w g z \right)\f$ is the Darcy velocity,
 * - \f$D\f$ is the hydrodynamic dispersion tensor,
 * - \f$\vartheta\f$ is the decay rate.
 *
 * For the hydrodynamic dispersion tensor the relation
 * \f[
 * D = (\phi D_d + \beta_T \|\vec{q}\|) I + (\beta_L - \beta_T) \frac{\vec{q}
 * \vec{q}^T}{\|\vec{q}\|}
 * \f]
 * is implemented, where \f$D_d\f$ is the molecular diffusion coefficient,
 * \f$\beta_L\f$ the longitudinal dispersivity of chemical species, and
 * \f$\beta_T\f$ the transverse dispersivity of chemical species.
 *
 * The implementation uses a monolithic approach, i.e., both processes
 * are assembled within one global system of equations.
 *
 * ## Process Coupling
 *
 * The advective term of the concentration equation is given by the Richards
 * flow process, i.e., the concentration distribution depends on
 * darcy velocity of the Richards flow process. On the other hand the
 * concentration dependencies of the viscosity and density in the groundwater
 * flow couples the unsaturated H process to the C process.
 *
 * \note At the moment there is not any coupling by source or sink terms, i.e.,
 * the coupling is implemented only by density changes due to concentration
 * changes in the buoyance term in the groundwater flow. This coupling schema is
 * referred to as the Boussinesq approximation.
 * */
class RichardsComponentTransportProcess final : public Process
{
public:
    RichardsComponentTransportProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        RichardsComponentTransportProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller,
        bool const use_monolithic_scheme);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

private:
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

    RichardsComponentTransportProcessData _process_data;

    std::vector<
        std::unique_ptr<RichardsComponentTransportLocalAssemblerInterface>>
        _local_assemblers;
};

}  // namespace RichardsComponentTransport
}  // namespace ProcessLib
