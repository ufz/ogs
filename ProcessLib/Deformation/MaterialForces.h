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
#include "ProcessLib/Deformation/GMatrix.h"
#include "ProcessLib/Parameter/SpatialPosition.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
namespace SmallDeformation
{
struct MaterialForcesInterface
{
    virtual std::vector<double> const& getMaterialForces(
        std::vector<double> const& local_x,
        std::vector<double>& nodal_values) = 0;

    virtual ~MaterialForcesInterface() = default;
};

template <int DisplacementDim, typename ShapeFunction,
          typename ShapeMatricesType, typename NodalForceVectorType,
          typename NodalDisplacementVectorType, typename GradientVectorType,
          typename GradientMatrixType, typename IPData,
          typename IntegrationMethod>
std::vector<double> const& getMaterialForces(
    std::vector<double> const& local_x, std::vector<double>& nodal_values,
    IntegrationMethod const& _integration_method, IPData const& _ip_data,
    MeshLib::Element const& element, bool const is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    nodal_values.clear();
    auto local_b = MathLib::createZeroedVector<NodalDisplacementVectorType>(
        nodal_values, ShapeFunction::NPOINTS * DisplacementDim);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sigma = _ip_data[ip].sigma;
        auto const& N = _ip_data[ip].N;
        auto const& dNdx = _ip_data[ip].dNdx;

        auto const& psi = _ip_data[ip].free_energy_density;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                element, _ip_data[ip].N);

        // For the 2D case the 33-component is needed (and the four entries
        // of the non-symmetric matrix); In 3d there are nine entries.
        GradientMatrixType G(
            DisplacementDim * DisplacementDim + (DisplacementDim == 2 ? 1 : 0),
            ShapeFunction::NPOINTS * DisplacementDim);
        Deformation::computeGMatrix<DisplacementDim, ShapeFunction::NPOINTS>(
            dNdx, G, is_axially_symmetric, N, x_coord);

        GradientVectorType const grad_u =
            G * Eigen::Map<NodalForceVectorType const>(
                    local_x.data(), ShapeFunction::NPOINTS * DisplacementDim);

        GradientVectorType eshelby_stress;
        eshelby_stress.setZero(DisplacementDim * DisplacementDim +
                               (DisplacementDim == 2 ? 1 : 0));

        if (DisplacementDim == 3)
        {
            eshelby_stress[0] = eshelby_stress[DisplacementDim + 1] =
                eshelby_stress[8] = psi;

            eshelby_stress[0] -= sigma[0] * grad_u[0] +
                                 sigma[3] / std::sqrt(2) * grad_u[3] +
                                 sigma[5] / std::sqrt(2) * grad_u[6];

            eshelby_stress[1] -= sigma[3] / std::sqrt(2) * grad_u[0] +
                                 sigma[1] * grad_u[3] +
                                 sigma[4] / std::sqrt(2) * grad_u[6];

            eshelby_stress[2] -= sigma[5] / std::sqrt(2) * grad_u[0] +
                                 sigma[4] / std::sqrt(2) * grad_u[3] +
                                 sigma[2] * grad_u[6];

            eshelby_stress[3] -= sigma[0] * grad_u[1] +
                                 sigma[3] / std::sqrt(2) * grad_u[4] +
                                 sigma[5] / std::sqrt(2) * grad_u[7];

            eshelby_stress[4] -= sigma[3] / std::sqrt(2) * grad_u[1] +
                                 sigma[1] * grad_u[4] +
                                 sigma[4] / std::sqrt(2) * grad_u[7];

            eshelby_stress[5] -= sigma[5] / std::sqrt(2) * grad_u[1] +
                                 sigma[4] / std::sqrt(2) * grad_u[4] +
                                 sigma[2] * grad_u[7];

            eshelby_stress[6] -= sigma[0] * grad_u[2] +
                                 sigma[3] / std::sqrt(2) * grad_u[5] +
                                 sigma[5] / std::sqrt(2) * grad_u[8];

            eshelby_stress[7] -= sigma[3] / std::sqrt(2) * grad_u[2] +
                                 sigma[1] * grad_u[5] +
                                 sigma[4] / std::sqrt(2) * grad_u[8];

            eshelby_stress[8] -= sigma[5] / std::sqrt(2) * grad_u[2] +
                                 sigma[4] / std::sqrt(2) * grad_u[5] +
                                 sigma[2] * grad_u[8];
        }
        else if (DisplacementDim == 2)
        {
            eshelby_stress[0] = eshelby_stress[DisplacementDim + 1] =
                eshelby_stress[4] = psi;

            eshelby_stress[0] -=
                sigma[0] * grad_u[0] + sigma[3] / std::sqrt(2) * grad_u[2];

            eshelby_stress[1] -=
                sigma[3] / std::sqrt(2) * grad_u[0] + sigma[1] * grad_u[2];

            eshelby_stress[2] -=
                sigma[0] * grad_u[1] + sigma[3] / std::sqrt(2) * grad_u[3];

            eshelby_stress[3] -=
                sigma[3] / std::sqrt(2) * grad_u[1] + sigma[1] * grad_u[3];

            // for axial-symmetric case the following is not zero in general
            eshelby_stress[4] -= sigma[2] * grad_u[4];
        }
        else
            OGS_FATAL(
                "Material forces not implemented for displacement dimension "
                "other than 2 and 3.");

        auto const& w = _ip_data[ip].integration_weight;
        local_b += G.transpose() * eshelby_stress * w;
    }

    return nodal_values;
}

template <typename LocalAssemblerInterface>
void writeMaterialForces(
    MeshLib::PropertyVector<double>& material_forces,
    std::vector<std::unique_ptr<LocalAssemblerInterface>> const&
        local_assemblers,
    NumLib::LocalToGlobalIndexMap const& local_to_global_index_map,
    GlobalVector const& x)
{
    DBUG("Compute material forces for small deformation process.");

    // Material forces
    std::fill(std::begin(material_forces), std::end(material_forces), 0);

    GlobalExecutor::executeDereferenced(
        [](const std::size_t mesh_item_id,
           LocalAssemblerInterface& local_assembler,
           const NumLib::LocalToGlobalIndexMap& dof_table,
           GlobalVector const& x,
           std::vector<double>& node_values) {
            auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
            std::vector<double> local_data;
            auto const local_x = x.get(indices);

            local_assembler.getMaterialForces(local_x, local_data);

            assert(local_data.size() == indices.size());
            for (std::size_t i = 0; i < indices.size(); ++i)
                node_values[indices[i]] += local_data[i];
        },
        local_assemblers, local_to_global_index_map, x, material_forces);
}

}  // namespace SmallDeformation
}  // namespace ProcessLib
