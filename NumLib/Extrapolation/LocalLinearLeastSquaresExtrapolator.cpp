/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LocalLinearLeastSquaresExtrapolator.h"

#include <functional>

#include <Eigen/Core>
#include <logog/include/logog.hpp>

#include "MathLib/LinAlg/LinAlg.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "NumLib/Assembler/SerialExecutor.h"
#include "NumLib/Function/Interpolation.h"
#include "ExtrapolatableElementCollection.h"

namespace NumLib
{
void LocalLinearLeastSquaresExtrapolator::extrapolate(
    ExtrapolatableElementCollection const& extrapolatables)
{
    _nodal_values.setZero();

    // counts the writes to each nodal value, i.e., the summands in order to
    // compute the average afterwards
    auto counts =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(_nodal_values);
    counts->setZero();  // TODO BLAS?

    auto const size = extrapolatables.size();
    for (std::size_t i=0; i<size; ++i) {
        extrapolateElement(i, extrapolatables, *counts);
    }

    MathLib::LinAlg::componentwiseDivide(_nodal_values, _nodal_values, *counts);
}

void LocalLinearLeastSquaresExtrapolator::calculateResiduals(
    ExtrapolatableElementCollection const& extrapolatables)
{
    assert(static_cast<std::size_t>(_residuals.size()) ==
           extrapolatables.size());

    auto const size = extrapolatables.size();
    for (std::size_t i=0; i<size; ++i) {
        calculateResiudalElement(i, extrapolatables);
    }
}

void LocalLinearLeastSquaresExtrapolator::extrapolateElement(
    std::size_t const element_index,
    ExtrapolatableElementCollection const& extrapolatables,
    GlobalVector& counts)
{
    auto const& integration_point_values =
        extrapolatables.getIntegrationPointValues(
            element_index, _integration_point_values_cache);

    // number of nodes in the element
    const auto nn = extrapolatables.getShapeMatrix(element_index, 0).cols();
    // number of integration points in the element
    const auto ni = integration_point_values.size();

    assert(ni >= static_cast<decltype(ni)>(nn) &&
           "Least squares is not possible if there are more nodes than"
           "integration points.");

    auto& N = _local_matrix_cache; // TODO make that local?
    N.resize(ni, nn); // TODO: might reallocate very often

    for (auto int_pt = decltype(ni){0}; int_pt < ni; ++int_pt) {
        auto const& shp_mat =
            extrapolatables.getShapeMatrix(element_index, int_pt);
        assert(shp_mat.cols() == nn);

        // copy shape matrix to extrapolation matrix columnwise
        N.block(int_pt, 0, 1, nn) = shp_mat;
    }

    // TODO make gp_vals an Eigen::VectorXd const& ?
    Eigen::Map<const Eigen::VectorXd> const integration_point_values_vec(
        integration_point_values.data(), integration_point_values.size());

    // TODO
    // optimization: Store decomposition of N*N^T or N^T for reuse?
    //   cf. http://eigen.tuxfamily.org/dox/classEigen_1_1FullPivLU.html
    //   cf. http://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html
    //   cf. http://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html
    //   cf. https://eigen.tuxfamily.org/dox/classEigen_1_1HouseholderQR.html
    // Extrapolate several values at once?

    // do the least squares computation using QR decomposition.
    Eigen::VectorXd tmp = N.householderQr().solve(integration_point_values_vec);

    // option: do the least squares computation using LDLT decomposition.
    // Eigen::VectorXd tmp =
    //         (N*N.transpose()).ldlt().solve(N*integration_point_values_vec);

    // TODO: for now always zeroth component is used
    auto const& global_indices = _local_to_global(element_index, 0).rows;

    _nodal_values.add(global_indices,
                      tmp);  // TODO does that give rise to PETSc problems?
    counts.add(global_indices, std::vector<double>(global_indices.size(), 1.0));
}

void LocalLinearLeastSquaresExtrapolator::calculateResiudalElement(
    std::size_t const element_index,
    ExtrapolatableElementCollection const& extrapolatables)
{
    auto const& gp_vals = extrapolatables.getIntegrationPointValues(
        element_index, _integration_point_values_cache);
    const unsigned ni = gp_vals.size();  // number of gauss points

    // TODO: for now always zeroth component is used
    const auto& global_indices = _local_to_global(element_index, 0).rows;

    // filter nodal values of the current element
    std::vector<double> nodal_vals_element;
    nodal_vals_element.resize(global_indices.size());
    for (unsigned i = 0; i < global_indices.size(); ++i) {
        // TODO PETSc negative indices?
        nodal_vals_element[i] = _nodal_values[global_indices[i]];
    }

    double residual = 0.0;
    double gp_val_extrapol = 0.0;

    for (unsigned gp = 0; gp < ni; ++gp) {
        NumLib::shapeFunctionInterpolate(
            nodal_vals_element,
            extrapolatables.getShapeMatrix(element_index, gp),
            gp_val_extrapol);
        auto const& ax_m_b = gp_val_extrapol - gp_vals[gp];
        residual += ax_m_b * ax_m_b;
    }

    _residuals.set(element_index, std::sqrt(residual / ni));
}

}  // namespace NumLib
