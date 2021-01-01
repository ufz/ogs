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

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Process.h"
#include "SteadyStateDiffusionFEM.h"
#include "SteadyStateDiffusionData.h"
#include "ProcessLib/SurfaceFlux/SurfaceFluxData.h"

namespace ProcessLib
{
namespace SteadyStateDiffusion
{
class SteadyStateDiffusion final : public Process
{
public:
    SteadyStateDiffusion(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        SteadyStateDiffusionData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return true; }
    //! @}

    Eigen::Vector3d getFlux(std::size_t element_id,
                            MathLib::Point3d const& p,
                            double const t,
                            std::vector<GlobalVector*> const& x) const override
    {
        // fetch local_x from primary variable
        std::vector<GlobalIndexType> indices_cache;
        auto const r_c_indices = NumLib::getRowColumnIndices(
            element_id, *_local_to_global_index_map, indices_cache);
        constexpr int process_id = 0;  // monolithic scheme.
        std::vector<double> local_x(x[process_id]->get(r_c_indices.rows));

        return _local_assemblers[element_id]->getFlux(p, t, local_x);
    }

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     const double t,
                                     const double /*delta_t*/,
                                     int const process_id) override
    {
        // For this single process, process_id is always zero.
        if (process_id != 0)
        {
            OGS_FATAL(
                "The condition of process_id = 0 must be satisfied for "
                "SteadyStateDiffusion, which is a single process.");
        }
        if (!_surfaceflux)  // computing the surfaceflux is optional
        {
            return;
        }

        ProcessLib::ProcessVariable const& pv =
            getProcessVariables(process_id)[0];

        _surfaceflux->integrate(x, t, *this, process_id, _integration_order,
                                _mesh, pv.getActiveElementIDs());
    }

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

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

    SteadyStateDiffusionData _process_data;

    std::vector<std::unique_ptr<SteadyStateDiffusionLocalAssemblerInterface>>
        _local_assemblers;

    std::unique_ptr<ProcessLib::SurfaceFluxData> _surfaceflux;
};

}  // namespace SteadyStateDiffusion
}   // namespace ProcessLib
