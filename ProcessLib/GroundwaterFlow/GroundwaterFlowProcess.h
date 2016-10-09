/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
#define PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "ProcessLib/Process.h"
#include "GroundwaterFlowFEM.h"
#include "GroundwaterFlowProcessData.h"
#include "ProcessLib/CalculateSurfaceFlux/CalculateSurfaceFlux.h"

// TODO used for output, if output classes are ready this has to be changed
#include "MeshLib/IO/writeMeshToFile.h"

namespace ProcessLib
{
namespace GroundwaterFlow
{
class GroundwaterFlowProcess final : public Process
{
    using Base = Process;

public:
    GroundwaterFlowProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        GroundwaterFlowProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller,
        MeshLib::Mesh* balance_mesh, std::string&& balance_pv_name,
        std::string&& balance_out_fname);

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return true; }
    //! @}

    std::vector<double> getFlux(std::size_t element_id,
                                MathLib::Point3d const& p,
                                GlobalVector const& x) const override
    {
        // fetch local_x from primary variable
        std::vector<GlobalIndexType> indices_cache;
        auto const r_c_indices = NumLib::getRowColumnIndices(
            element_id, *_local_to_global_index_map, indices_cache);
        std::vector<double> local_x(x.get(r_c_indices.rows));

        return _local_assemblers[element_id]->getFlux(p, local_x);
    }

    void postTimestepConcreteProcess(GlobalVector const& x) override
    {
        if (_balance_mesh) // computing the balance is optional
        {
            std::vector<double> init_values(
                _balance_mesh->getNumberOfElements(), 0.0);
            MeshLib::addPropertyToMesh(*_balance_mesh, _balance_pv_name,
                                       MeshLib::MeshItemType::Cell, 1,
                                       init_values);
            auto balance = ProcessLib::CalculateSurfaceFlux(
                *_balance_mesh,
                getProcessVariables()[0].get().getNumberOfComponents(),
                _integration_order);

            auto* const balance_pv =
                _balance_mesh->getProperties()
                    .template getPropertyVector<double>(_balance_pv_name);

            balance.integrate(x, *balance_pv, *this);
            // post: surface_mesh has vectorial element property

            // TODO output, if output classes are ready this has to be
            // changed
            MeshLib::IO::writeMeshToFile(*_balance_mesh, _balance_out_fname);
        }
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

    std::unique_ptr<MeshLib::Mesh> _balance_mesh;
    std::string const _balance_pv_name;
    std::string const _balance_out_fname;
};

}   // namespace GroundwaterFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
