/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>

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
struct HTMaterialProperties;

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
 * \note At the moment there is not any coupling by source or sink terms, i.e.,
 * the coupling is implemented only by density changes due to temperature
 * changes in the buoyance term in the groundwater flow. This coupling schema is
 * referred to as the Boussinesq approximation.
 * */
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
        std::unique_ptr<HTMaterialProperties>&& material_properties,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller,
        bool const use_monolithic_scheme,
        std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux,
        const int heat_transport_process_id,
        const int hydraulic_process_id);
    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

    Eigen::Vector3d getFlux(std::size_t element_id,
                            MathLib::Point3d const& p,
                            double const t,
                            GlobalVector const& x) const override;

    void setCoupledTermForTheStaggeredSchemeToLocalAssemblers() override;

    void postTimestepConcreteProcess(GlobalVector const& x,
                                     const double t,
                                     const double delta_t,
                                     int const process_id) override;

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

    void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                                    double const dt,
                                    const int process_id) override;

    void setCoupledSolutionsOfPreviousTimeStepPerProcess(const int process_id);

    /// Set the solutions of the previous time step to the coupled term.
    /// It is only for the staggered scheme, and it must be called within
    /// the coupling loop because that the coupling term is only created there.
    void setCoupledSolutionsOfPreviousTimeStep();

    /**
     * @copydoc ProcessLib::Process::getDOFTableForExtrapolatorData()
     */
    std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
        getDOFTableForExtrapolatorData() const override;

    const std::unique_ptr<HTMaterialProperties> _material_properties;

    std::vector<std::unique_ptr<HTLocalAssemblerInterface>> _local_assemblers;

    /// Solutions of the previous time step
    std::array<std::unique_ptr<GlobalVector>, 2> _xs_previous_timestep;

    std::unique_ptr<ProcessLib::SurfaceFluxData> _surfaceflux;

    const int _heat_transport_process_id;
    const int _hydraulic_process_id;
};

}  // namespace HT
}  // namespace ProcessLib
