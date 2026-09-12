// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <array>

#include "HTProcessData.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "ProcessLib/Process.h"

namespace NumLib
{
class LocalToGlobalIndexMap;
}

namespace ProcessLib
{
struct SurfaceFluxData;

namespace HT
{
class HTLocalAssemblerInterface;

/**
 * # HT process
 *
 * The implementation uses a monolithic approach, i.e., both processes
 * are assembled within one global system of equations.
 *
 * ## Process Coupling
 *
 * The advective term of the heat conduction equation is given by the confined
 * groundwater flow process, i.e., the heat conduction depends on darcy velocity
 * of the groundwater flow process. On the other hand the temperature
 * dependencies of the viscosity and density in the groundwater flow couples the
 * H process to the T process.
 *
 * \note
 *    - At the moment there is no coupling by source or sink terms, i.e., the
 *      coupling is implemented only through density changes due to
 *      temperature changes in the buoyancy term of the groundwater flow.
 *      This coupling scheme is referred to as the Boussinesq approximation.
 *    - The fluid phase contribution to the storage coefficient is computed
 *      from the fluid compressibility, i.e.
 *      \f$\phi\frac{\partial \varrho_f}{\partial p}/\varrho_f\f$ with
 *      \f$\phi\f$ the porosity, \f$p\f$ the pore pressure, and
 *      \f$\varrho_f\f$ the fluid density.
 *    - The storage input parameter is for the solid phase only, and can be
 *      computed from the Biot coefficient \f$\alpha_B\f$ and the drained bulk
 *      modulus \f$K\f$ as \f$(\alpha_B-\phi)(1-\alpha_B)/K\f$. Equivalently,
 *      it is \f$(\alpha_B-\phi)/K_s\f$ with
 *      \f$K_s=K/(1-\alpha_B)\f$ the intrinsic bulk modulus of the solid
 *      phase. Therefore, if the Biot coefficient is defined as one, the
 *      storage input parameter must be zero. This is enforced whenever
 *      \f$\alpha_T^s\f$ is defined, since that is when \f$\alpha_B\f$ is read
 *      at all: checkBiotStorageRelation() compares the evaluated values at
 *      the integration points of every element once during initialisation, at
 *      \f$t=0\f$, which covers all properties independent of the primary
 *      variables and of the time, and
 *      evalEffectiveThermalExpansivity() compares them again at each
 *      integration point during assembly, which covers the remaining property
 *      types.
 *    - The input parameters of the Biot coefficient \f$\alpha_B\f$ and the
 *      solid thermal expansivity (linear) \f$\alpha_T^s\f$ are optional. Only
 *      one direction is enforced: if \f$\alpha_T^s\f$ is given, then
 *      \f$\alpha_B\f$ must be given too, see
 *      checkThermalExpansivitySetting(). The reverse is not enforced, because
 *      \f$\alpha_B\f$ is never read without \f$\alpha_T^s\f$; it is then
 *      silently ignored. They are only used to compute the effective thermal
 *      expansivity, which is defined as:
 *      \f[
 *          3(\alpha_B-\phi)\alpha_T^s - \phi \frac{\partial \varrho_f}
 *          {\partial T}/\varrho_f
 *      \f]
 *      If they are not defined, the effective thermal expansivity is computed
 *      as \f$-\phi \frac{\partial \varrho_f}{\partial T}/\varrho_f\f$.
 *    - The storage term of the pressure equation is
 *      \f$\phi\frac{\partial \varrho_f}{\partial p}/\varrho_f + S_s\f$ with
 *      \f$S_s\f$ the storage input parameter. It does not contain the thermal
 *      expansivity, and therefore it vanishes exactly when the liquid density
 *      does not depend on the pressure and \f$S_s\f$ is zero -- for instance
 *      for a temperature-only density model combined with \f$\alpha_B=1\f$,
 *      which forces \f$S_s=0\f$. Such a setup is physically inconsistent and
 *      numerical instability can occur.
 *    - The governing equation can be set to either a volume balance or a mass
 *      balance. The default is a volume balance. If the governing equation is
 *      set to a mass balance, the input of the fluid phase boundary and
 *      source/sink terms changes from a volume rate to a mass rate: the unit
 *      of the Neumann boundary condition is changed from \f$[m/s]\f$ to
 *      \f$[kg/(m^2\,s)]\f$, and the unit of the source/sink term from
 *      \f$[m^3/s]\f$ to \f$[kg/s]\f$.
 */
class HTProcess final : public Process
{
public:
    HTProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        HTProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        bool const use_monolithic_scheme,
        std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux);
    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

    Eigen::Vector3d getFlux(std::size_t element_id,
                            MathLib::Point3d const& p,
                            double const t,
                            std::vector<GlobalVector*> const& x) const override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     std::vector<GlobalVector*> const& x_prev,
                                     const double t,
                                     const double delta_t,
                                     int const process_id) override;

private:
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

    /**
     * @copydoc ProcessLib::Process::getDOFTableForExtrapolatorData()
     */
    std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
    getDOFTableForExtrapolatorData() const override;

    HTProcessData _process_data;

    std::vector<std::unique_ptr<HTLocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<ProcessLib::SurfaceFluxData> _surfaceflux;
};

}  // namespace HT
}  // namespace ProcessLib
