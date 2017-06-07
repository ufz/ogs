/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>
#include <vector>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
namespace SmallDeformation
{
struct NodalForceCalculationInterface
{
    virtual std::vector<double> const& getNodalForces(
        std::vector<double>& nodal_values) const = 0;

    virtual ~NodalForceCalculationInterface() = default;
};

template <int DisplacementDim, typename ShapeFunction,
          typename ShapeMatricesType, typename NodalDisplacementVectorType,
          typename BMatrixType, typename IPData, typename IntegrationMethod>
std::vector<double> const& getNodalForces(
    std::vector<double>& nodal_values,
    IntegrationMethod const& _integration_method, IPData const& _ip_data,
    MeshLib::Element const& element, bool const is_axially_symmetric)
{
    nodal_values.clear();
    auto local_b = MathLib::createZeroedVector<NodalDisplacementVectorType>(
        nodal_values, ShapeFunction::NPOINTS * DisplacementDim);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition x_position;
    x_position.setElementID(element.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        BMatrixType B(KelvinVectorDimensions<DisplacementDim>(),
                      ShapeFunction::NPOINTS * DisplacementDim);
        auto const x_coord =
            interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                element, _ip_data[ip].N);
        LinearBMatrix::computeBMatrix<DisplacementDim, ShapeFunction::NPOINTS>(
            _ip_data[ip].dNdx, B, is_axially_symmetric, _ip_data[ip].N,
            x_coord);
        auto& sigma = _ip_data[ip].sigma;

        local_b.noalias() += B.transpose() * sigma * w;
    }

    return nodal_values;
}

template <typename LocalAssemblerInterface>
void writeNodalForces(
    MeshLib::PropertyVector<double>& nodal_forces,
    std::vector<std::unique_ptr<LocalAssemblerInterface>> const&
        local_assemblers,
    NumLib::LocalToGlobalIndexMap const& local_to_global_index_map)
{
    DBUG("Compute nodal forces for small deformation process.");

    // Zero-out the output vector before averaging.
    std::fill(std::begin(nodal_forces), std::end(nodal_forces), 0);

    GlobalExecutor::executeDereferenced(
        [](const std::size_t mesh_item_id,
           LocalAssemblerInterface& local_assembler,
           const NumLib::LocalToGlobalIndexMap& dof_table,
           std::vector<double>& node_values) {
            auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
            std::vector<double> local_data;

            local_assembler.getNodalForces(local_data);

            assert(local_data.size() == indices.size());
            for (std::size_t i = 0; i < indices.size(); ++i)
                node_values[indices[i]] += local_data[i];
        },
        local_assemblers, local_to_global_index_map, nodal_forces);
}

}  // namespace SmallDeformation
}  // namespace ProcessLib
