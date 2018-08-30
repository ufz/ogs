/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "GroundwaterFlowFEM.h"
#include "GroundwaterFlowProcessData.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/CalculateSurfaceFlux/Balance.h"
#include "ProcessLib/Process.h"

// TODO used for output, if output classes are ready this has to be changed
#include "MeshLib/IO/writeMeshToFile.h"

namespace ProcessLib
{
namespace GroundwaterFlow
{
class GroundwaterFlowProcess final : public Process
{
public:
    GroundwaterFlowProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        GroundwaterFlowProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller,
        std::unique_ptr<ProcessLib::Balance>&& balance);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return true; }
    //! @}

    Eigen::Vector3d getFlux(std::size_t element_id,
                            MathLib::Point3d const& p,
                            double const t,
                            GlobalVector const& x) const override
    {
        // fetch local_x from primary variable
        std::vector<GlobalIndexType> indices_cache;
        auto const r_c_indices = NumLib::getRowColumnIndices(
            element_id, *_local_to_global_index_map, indices_cache);
        std::vector<double> local_x(x.get(r_c_indices.rows));

        return _local_assemblers[element_id]->getFlux(p, t, local_x);
    }

    void postTimestepConcreteProcess(GlobalVector const& x,
                                     const double t,
                                     const double /*delta_t*/,
                                     int const process_id) override
    {
        //For this single process, process_id is always zero.
        if (process_id != 0)
        {
            OGS_FATAL("The condition of process_id = 0 must be satisfied for "
                      "GroundwaterFlowProcess, which is a single process." );
        }
        if (!_balance)  // computing the balance is optional
        {
            return;
        }
        _balance->integrate(x, t, *this, process_id, _integration_order, _mesh);
        _balance->save(t);
    }

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

    GroundwaterFlowProcessData _process_data;

    std::vector<std::unique_ptr<GroundwaterFlowLocalAssemblerInterface>>
        _local_assemblers;

    std::unique_ptr<ProcessLib::Balance> _balance;
};

}   // namespace GroundwaterFlow
}   // namespace ProcessLib
