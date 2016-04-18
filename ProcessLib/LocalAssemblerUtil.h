/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_LOCALASSEMBLERUTIL_H
#define PROCESSLIB_LOCALASSEMBLERUTIL_H

#include<vector>
#include<Eigen/Core>

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"

namespace MeshLib { class Element; }

namespace ProcessLib
{

template<typename ShapeFunction, typename ShapeMatricesType, typename IntegrationMethod,
         unsigned GlobalDim>
std::vector<typename ShapeMatricesType::ShapeMatrices>
initShapeMatrices(MeshLib::Element const& e, unsigned integration_order)
{
    std::vector<typename ShapeMatricesType::ShapeMatrices> shape_matrices;

    using FemType = NumLib::TemplateIsoparametric<
        ShapeFunction, ShapeMatricesType>;

    FemType fe(*static_cast<const typename ShapeFunction::MeshElement*>(&e));

    IntegrationMethod integration_method(integration_order);
    std::size_t const n_integration_points = integration_method.getNPoints();

    shape_matrices.reserve(n_integration_points);
    for (std::size_t ip = 0; ip < n_integration_points; ++ip) {
        shape_matrices.emplace_back(ShapeFunction::DIM, GlobalDim,
                                     ShapeFunction::NPOINTS);
        fe.computeShapeFunctions(
                integration_method.getWeightedPoint(ip).getCoords(),
                shape_matrices[ip]);
    }

    return shape_matrices;
}

Eigen::Map<Eigen::MatrixXd>
setupLocalMatrix(std::vector<double>& matrix_data, std::size_t const local_matrix_size);

Eigen::Map<Eigen::VectorXd>
setupLocalVector(std::vector<double>& matrix_data, std::size_t const local_matrix_size);

} // namespace ProcessLib

#endif // PROCESSLIB_LOCALASSEMBLERUTIL_H
