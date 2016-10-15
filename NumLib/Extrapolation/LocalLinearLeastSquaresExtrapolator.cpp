/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LocalLinearLeastSquaresExtrapolator.h"

#include <Eigen/SVD>
#include <logog/include/logog.hpp>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "NumLib/Assembler/SerialExecutor.h"
#include "NumLib/Function/Interpolation.h"
#include "ExtrapolatableElementCollection.h"

namespace NumLib
{
LocalLinearLeastSquaresExtrapolator::LocalLinearLeastSquaresExtrapolator(
    NumLib::LocalToGlobalIndexMap const& dof_table)
    : _nodal_values(NumLib::GlobalVectorProvider::provider.getVector(
          MathLib::MatrixSpecifications(dof_table.dofSizeWithoutGhosts(),
                                        dof_table.dofSizeWithoutGhosts(),
                                        &dof_table.getGhostIndices(),
                                        nullptr)))
#ifndef USE_PETSC
    , _residuals(dof_table.size())
#else
    , _residuals(dof_table.size(), false)
#endif
    , _local_to_global(dof_table)
{
    /* Note in case the following assertion fails:
     * If you copied the extrapolation code, for your processes from
     * somewhere, note that the code from the groundwater flow process might
     * not suit your needs: It is a special case and is therefore most
     * likely too simplistic. You better adapt the extrapolation code from
     * some more advanced process, like the TES process.
     */
    assert(dof_table.getNumberOfComponents() == 1 &&
           "The d.o.f. table passed must be for one variable that has "
           "only one component!");
}

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
        calculateResidualElement(i, extrapolatables);
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

    auto const& N_0 = extrapolatables.getShapeMatrix(element_index, 0);
    const unsigned num_nodes = N_0.cols();
    const unsigned num_int_pts = integration_point_values.size();

    assert(num_int_pts >= num_nodes &&
           "Least squares is not possible if there are more nodes than"
           "integration points.");

    auto const pair_it_inserted = _qr_decomposition_cache.emplace(
        std::make_pair(num_nodes, num_int_pts), CachedData{N_0});

    auto& cached_data = pair_it_inserted.first->second;
    if (pair_it_inserted.second)
    {
        DBUG("Computing new singular value decomposition");

        // interpolation_matrix * nodal_values = integration_point_values
        // We are going to pseudo-invert this relation now using singular value
        // decomposition.
        Eigen::MatrixXd interpolation_matrix(num_int_pts, num_nodes);

        for (unsigned int_pt = 0; int_pt < num_int_pts; ++int_pt) {
            auto const& shp_mat =
                extrapolatables.getShapeMatrix(element_index, int_pt);
            assert(shp_mat.cols() == num_nodes);

            // copy shape matrix to extrapolation matrix row-wise
            interpolation_matrix.row(int_pt) = shp_mat;
        }

        // JacobiSVD is extremely reliable, but fast only for small matrices.
        // But we usually have small matrices and we don't compute very often.
        // Cf. http://eigen.tuxfamily.org/dox/group__TopicLinearAlgebraDecompositions.html
        //
        // Decomposes interpolation_matrix = U S V^T.
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(
            interpolation_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);

        auto const& S = svd.singularValues();
        auto const& U = svd.matrixU();
        auto const& V = svd.matrixV();

        // Compute and save the pseudo inverse V * S^{-1} * U^T.
        auto const rank = svd.rank();
        assert(rank == num_nodes);

        // cf. http://eigen.tuxfamily.org/dox/JacobiSVD_8h_source.html
        cached_data.p_inv.noalias() =
            V.leftCols(rank) *
            S.head(rank).asDiagonal().inverse() *
            U.leftCols(rank).transpose();
    }
    else if (cached_data.N_0 != N_0) {
        OGS_FATAL("The cached and the passed shapematrices differ.");
    }

    // TODO make gp_vals an Eigen::VectorXd const& ?
    auto const integration_point_values_vec =
        MathLib::toVector(integration_point_values);

    // Apply the pre-computed pseudo-inverse.
    Eigen::VectorXd const nodal_values =
        cached_data.p_inv * integration_point_values_vec;

    // TODO: for now always zeroth component is used. This has to be extended if
    // multi-component properties shall be extrapolated
    auto const& global_indices = _local_to_global(element_index, 0).rows;

    // TODO does that give rise to PETSc problems?
    _nodal_values.add(global_indices, nodal_values);
    counts.add(global_indices, std::vector<double>(global_indices.size(), 1.0));
}

void LocalLinearLeastSquaresExtrapolator::calculateResidualElement(
    std::size_t const element_index,
    ExtrapolatableElementCollection const& extrapolatables)
{
    auto const& int_pt_vals = extrapolatables.getIntegrationPointValues(
        element_index, _integration_point_values_cache);

    // TODO: for now always zeroth component is used
    const auto& global_indices = _local_to_global(element_index, 0).rows;

    const unsigned num_int_pts = int_pt_vals.size();
    const unsigned num_nodes = global_indices.size();

    // filter nodal values of the current element
    std::vector<double> nodal_vals_element(num_nodes);
    for (unsigned i = 0; i < num_nodes; ++i) {
        // TODO PETSc negative indices?
        nodal_vals_element[i] = _nodal_values[global_indices[i]];
    }

    double residual = 0.0;
    double int_pt_val_extrapolated = 0.0;

    for (unsigned int_pt = 0; int_pt < num_int_pts; ++int_pt) {
        NumLib::shapeFunctionInterpolate(
            nodal_vals_element,
            extrapolatables.getShapeMatrix(element_index, int_pt),
            int_pt_val_extrapolated);
        auto const& ax_m_b = int_pt_val_extrapolated - int_pt_vals[int_pt];
        residual += ax_m_b * ax_m_b;
    }

    _residuals.set(element_index, std::sqrt(residual / num_int_pts));
}

}  // namespace NumLib
