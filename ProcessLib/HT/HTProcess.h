/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
        HTProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
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
                            std::vector<GlobalVector*> const& x) const override;

    void setCoupledTermForTheStaggeredSchemeToLocalAssemblers(
        int const process_id) override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     const double t,
                                     const double delta_t,
                                     int const process_id) override;

private:
    template <unsigned GlobalDim>
    void checkProperties(MeshLib::Mesh const& mesh) const
    {
        // only needed as dummy for checking of existence of properties
        MaterialPropertyLib::VariableArray vars;
        double const t = 0.0;

        DBUG("Check the media properties ...");
        for (auto const& element : mesh.getElements())
        {
            auto const element_id = element->getID();
            ParameterLib::SpatialPosition pos;
            pos.setElementID(element_id);

            // check if a definition of the porous media exists
            auto const& medium =
                *_process_data.media_map->getMedium(element_id);

            // checking general medium properties
            auto const porosity =
                medium.property(MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(vars, pos, t);
            auto const K = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium.property(MaterialPropertyLib::PropertyType::permeability)
                    .value(vars, pos, t));

            // check if liquid phase definition and the corresponding properties
            // exists
            auto const& liquid_phase = medium.phase("AqueousLiquid");
            auto const mu =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::viscosity)
                    .template value<double>(vars, pos, t);
            auto const liquid_density =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, t);
            auto const specific_heat_capacity_fluid =
                liquid_phase
                    .property(MaterialPropertyLib::specific_heat_capacity)
                    .template value<double>(vars, pos, t);

            // check if solid phase definition and the corresponding properties
            // exists
            auto const& solid_phase = medium.phase("Solid");
            auto const specific_heat_capacity_solid =
                solid_phase
                    .property(MaterialPropertyLib::PropertyType::
                                  specific_heat_capacity)
                    .template value<double>(vars, pos, t);
            auto const solid_density =
                solid_phase.property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, t);
            auto const specific_storage =
                solid_phase.property(MaterialPropertyLib::PropertyType::storage)
                    .template value<double>(vars, pos, t);
        }
        DBUG("Media properties verified.");
    }

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                    double const t, double const dt,
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

    HTProcessData _process_data;

    std::vector<std::unique_ptr<HTLocalAssemblerInterface>> _local_assemblers;

    /// Solutions of the previous time step
    std::array<std::unique_ptr<GlobalVector>, 2> _xs_previous_timestep;

    std::unique_ptr<ProcessLib::SurfaceFluxData> _surfaceflux;

    const int _heat_transport_process_id;
    const int _hydraulic_process_id;
};

}  // namespace HT
}  // namespace ProcessLib
