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

#include "ComponentTransportFEM.h"
#include "ComponentTransportProcessData.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "ProcessLib/Process.h"

namespace ChemistryLib
{
class ChemicalSolverInterface;
}

namespace ProcessLib
{
struct SurfaceFluxData;

namespace ComponentTransport
{
/**
 * # ComponentTransport process
 *
 * ## Governing equations
 *
 * The flow process is described by
 * \f[
 * \phi \frac{\partial \rho}{\partial p} \frac{\partial p}{\partial t}
 *     + \phi \frac{\partial \rho}{\partial C} \frac{\partial C}{\partial t}
 *     - \nabla \cdot \left[\frac{\kappa}{\mu(C)} \rho \nabla \left( p + \rho g
 * z \right)\right]
 *     + Q_p = 0,
 * \f]
 * where the storage \f$S\f$ has been substituted by
 *      \f$\phi \frac{\partial \rho}{\partial p}\f$,
 * \f$\phi\f$ is the porosity,
 * \f$C\f$ is the concentration,
 * \f$p\f$ is the pressure,
 * \f$\kappa\f$ is permeability,
 * \f$\mu\f$ is viscosity of the fluid,
 * \f$\rho\f$ is the density of the fluid, and
 * \f$g\f$ is the gravitational acceleration.
 *
 * The mass transport process is described by
 * \f[
 * \phi R C \frac{\partial \rho}{\partial p} \frac{\partial p}{\partial t}
 *     + \phi R \left(\rho + C \frac{\partial \rho}{\partial C}\right)
 * \frac{\partial C}{\partial t}
 *     - \nabla \cdot \left[\frac{\kappa}{\mu(C)} \rho C \nabla \left( p + \rho
 * g z \right)
 *          + \rho D \nabla C\right]
 *     + Q_C + R \vartheta \phi \rho C = 0,
 *
 * \f]
 * where \f$R\f$ is the retardation factor,
 * \f$\vec{q} = -\frac{\kappa}{\mu(C)} \nabla \left( p + \rho g z \right)\f$
 *      is the Darcy velocity,
 * \f$D\f$ is the hydrodynamic dispersion tensor,
 * \f$\vartheta\f$ is the decay rate.
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
 * The advective term of the concentration equation is given by the confined
 * groundwater flow process, i.e., the concentration distribution depends on
 * Darcy velocity of the groundwater flow process. On the other hand the
 * concentration dependencies of the viscosity and density in the groundwater
 * flow couples the H process to the C process.
 *
 * \note At the moment there is not any coupling by source or sink terms, i.e.,
 * the coupling is implemented only by density and viscosity changes due to
 * concentration changes as well as by the temporal derivatives of each
 * variable.
 * */
class ComponentTransportProcess final : public Process
{
public:
    ComponentTransportProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        ComponentTransportProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        bool const use_monolithic_scheme,
        std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux,
        std::unique_ptr<ChemistryLib::ChemicalSolverInterface>&&
            chemical_solver_interface);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

    Eigen::Vector3d getFlux(std::size_t const element_id,
                            MathLib::Point3d const& p, double const t,
                            std::vector<GlobalVector*> const& x) const override;

    void setCoupledTermForTheStaggeredSchemeToLocalAssemblers(
        int const process_id) override;

    void solveReactionEquation(std::vector<GlobalVector*>& x, double const t,
                               double const dt) override;

    std::vector<GlobalVector> interpolateNodalValuesToIntegrationPoints(
        std::vector<GlobalVector*> const& nodal_values_vectors) const override;

    void extrapolateIntegrationPointValuesToNodes(
        const double t,
        std::vector<GlobalVector*> const& integration_point_values_vectors,
        std::vector<GlobalVector*>& nodal_values_vectors) override;

    void computeSecondaryVariableConcrete(double const /*t*/,
                                          double const /*dt*/,
                                          std::vector<GlobalVector*> const& x,
                                          GlobalVector const& /*x_dot*/,
                                          int const /*process_id*/) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     const double t,
                                     const double dt,
                                     int const process_id) override;

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

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

    ComponentTransportProcessData _process_data;

    std::vector<std::unique_ptr<ComponentTransportLocalAssemblerInterface>>
        _local_assemblers;

    std::unique_ptr<ProcessLib::SurfaceFluxData> _surfaceflux;

    std::unique_ptr<ChemistryLib::ChemicalSolverInterface>
        _chemical_solver_interface;
};

}  // namespace ComponentTransport
}  // namespace ProcessLib
