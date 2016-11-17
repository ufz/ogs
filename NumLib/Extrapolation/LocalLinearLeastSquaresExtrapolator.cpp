/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
    : _dof_table_single_component(dof_table)
{
    /* Note in case the following assertion fails:
     * If you copied the extrapolation code, for your processes from
     * somewhere, note that the code from the groundwater flow process might
     * not suit your needs: It is a special case and is therefore most
     * likely too simplistic. You better adapt the extrapolation code from
     * some more advanced process, like the TES process.
     */
    if (dof_table.getNumberOfComponents() != 1)
        OGS_FATAL(
            "The d.o.f. table passed must be for one variable that has "
            "only one component!");
}

void LocalLinearLeastSquaresExtrapolator::extrapolate(
    const unsigned num_components,
    ExtrapolatableElementCollection const& extrapolatables,
    const double t,
    GlobalVector const& current_solution,
    LocalToGlobalIndexMap const& dof_table)
{
    auto const num_nodal_dof_result =
        _dof_table_single_component.dofSizeWithoutGhosts() * num_components;
    if (!_nodal_values ||
        static_cast<std::size_t>(_nodal_values->size()) != num_nodal_dof_result)
    {
        _nodal_values = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
            MathLib::MatrixSpecifications{
                num_nodal_dof_result, num_nodal_dof_result, nullptr, nullptr});
    }
    _nodal_values->setZero();

    // counts the writes to each nodal value, i.e., the summands in order to
    // compute the average afterwards
    auto counts =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(*_nodal_values);
    counts->setZero();  // TODO BLAS?

    auto const size = extrapolatables.size();
    for (std::size_t i = 0; i < size; ++i)
    {
        extrapolateElement(i, num_components, extrapolatables, t,
                           current_solution, dof_table, *counts);
    }

    MathLib::LinAlg::componentwiseDivide(*_nodal_values, *_nodal_values,
                                         *counts);
}

void LocalLinearLeastSquaresExtrapolator::calculateResiduals(
    const unsigned num_components,
    ExtrapolatableElementCollection const& extrapolatables,
    const double t,
    GlobalVector const& current_solution,
    LocalToGlobalIndexMap const& dof_table)
{
    auto const num_element_dof_result = static_cast<GlobalIndexType>(
        _dof_table_single_component.size() * num_components);

    if (!_residuals || _residuals->size() != num_element_dof_result)
    {
#ifndef USE_PETSC
        _residuals.reset(new GlobalVector{num_element_dof_result});
#else
        _residuals.reset(new GlobalVector{num_element_dof_result, false});
#endif
    }

    if (static_cast<std::size_t>(num_element_dof_result) !=
        extrapolatables.size() * num_components)
    {
        OGS_FATAL("mismatch in number of D.o.F.");
    }

    auto const size = extrapolatables.size();
    for (std::size_t i = 0; i < size; ++i)
    {
        calculateResidualElement(i, num_components, extrapolatables, t,
                                 current_solution, dof_table);
    }
}

void LocalLinearLeastSquaresExtrapolator::extrapolateElement(
    std::size_t const element_index,
    const unsigned num_components,
    ExtrapolatableElementCollection const& extrapolatables,
    const double t,
    GlobalVector const& current_solution,
    LocalToGlobalIndexMap const& dof_table,
    GlobalVector& counts)
{
    auto const& integration_point_values =
        extrapolatables.getIntegrationPointValues(
            element_index, t, current_solution, dof_table,
            _integration_point_values_cache);

    auto const& N_0 = extrapolatables.getShapeMatrix(element_index, 0);
    const unsigned num_nodes = N_0.cols();
    const unsigned num_values = integration_point_values.size();

    if (num_values % num_components != 0)
        OGS_FATAL(
            "The number of computed integration point values is not divisable "
            "by the number of num_components. Maybe the computed property is "
            "not a %d-component vector for each integration point.",
            num_components);

    // number of integration points in the element
    const auto num_int_pts = num_values / num_components;

    if (num_int_pts < static_cast<decltype(num_int_pts)>(num_nodes))
        OGS_FATAL(
            "Least squares is not possible if there are more nodes than"
            "integration points.");

    auto const pair_it_inserted = _qr_decomposition_cache.emplace(
        std::make_pair(num_nodes, num_int_pts), CachedData{});

    auto& cached_data = pair_it_inserted.first->second;
    if (pair_it_inserted.second)
    {
        DBUG("Computing new singular value decomposition");

        // interpolation_matrix * nodal_values = integration_point_values
        // We are going to pseudo-invert this relation now using singular value
        // decomposition.
        auto& interpolation_matrix = cached_data.A;
        interpolation_matrix.resize(num_int_pts, num_nodes);

        interpolation_matrix.row(0) = N_0;
        for (unsigned int_pt = 1; int_pt < num_int_pts; ++int_pt) {
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
        cached_data.A_pinv.noalias() =
            V.leftCols(rank) *
            S.head(rank).asDiagonal().inverse() *
            U.leftCols(rank).transpose();
    }
    else if (cached_data.A.row(0) != N_0) {
        OGS_FATAL("The cached and the passed shapematrices differ.");
    }

    // TODO: for now always zeroth component is used. This has to be extended if
    // multi-component properties shall be extrapolated
    auto const& global_indices =
        _dof_table_single_component(element_index, 0).rows;

    if (num_components == 1)
    {
        // TODO make gp_vals an Eigen::VectorXd const& ?
        auto const integration_point_values_vec =
            MathLib::toVector(integration_point_values);

        // Apply the pre-computed pseudo-inverse.
        Eigen::VectorXd const nodal_values =
            cached_data.A_pinv * integration_point_values_vec;

        // TODO does that give rise to PETSc problems?
        _nodal_values->add(global_indices, nodal_values);
        counts.add(global_indices,
                   std::vector<double>(global_indices.size(), 1.0));
    }
    else
    {
        // TODO make gp_vals an Eigen::VectorXd const& ?
        auto const integration_point_values_mat = MathLib::toMatrix(
            integration_point_values, num_components, num_int_pts);

        // Apply the pre-computed pseudo-inverse.
        Eigen::MatrixXd const nodal_values =
            cached_data.A_pinv * integration_point_values_mat.transpose();

        std::vector<GlobalIndexType> indices;
        indices.reserve(num_components * global_indices.size());

        // component-wise in nodal_values_flat, location-wise ordering in
        // _nodal_values
        for (unsigned comp = 0; comp < num_components; ++comp)
        {
            for (auto i : global_indices)
            {
                indices.push_back(num_components * i + comp);
            }
        }

        // Nodal_values are passed as a raw pointer, because PETScVector and
        // EigenVector implementations differ slightly.
        _nodal_values->add(indices, nodal_values.data());
        counts.add(indices, std::vector<double>(indices.size(), 1.0));
    }
}

void LocalLinearLeastSquaresExtrapolator::calculateResidualElement(
    std::size_t const element_index,
    const unsigned num_components,
    ExtrapolatableElementCollection const& extrapolatables,
    const double t,
    GlobalVector const& current_solution,
    LocalToGlobalIndexMap const& dof_table)
{
    auto const& int_pt_vals = extrapolatables.getIntegrationPointValues(
        element_index, t, current_solution, dof_table,
        _integration_point_values_cache);

    const unsigned num_values = int_pt_vals.size();
    if (num_values % num_components != 0)
        OGS_FATAL(
            "The number of computed integration point values is not divisable "
            "by the number of num_components. Maybe the computed property is "
            "not a %d-component vector for each integration point.",
            num_components);

    // number of integration points in the element
    const auto num_int_pts = num_values / num_components;

    // TODO: for now always zeroth component is used
    const auto& global_indices =
        _dof_table_single_component(element_index, 0).rows;
    const unsigned num_nodes = global_indices.size();

    auto const& interpolation_matrix =
        _qr_decomposition_cache.find({num_nodes, num_int_pts})->second.A;

    Eigen::VectorXd nodal_vals_element(num_nodes);
    auto const int_pt_vals_mat =
        MathLib::toMatrix(int_pt_vals, num_components, num_int_pts);

    for (unsigned comp = 0; comp < num_components; ++comp)
    {
        // filter nodal values of the current element
        for (unsigned i = 0; i < num_nodes; ++i)
        {
            // TODO PETSc negative indices?
            auto const idx = num_components * global_indices[i] + comp;
            nodal_vals_element[i] = (*_nodal_values)[idx];
        }

        double const residual = (interpolation_matrix * nodal_vals_element -
                                 int_pt_vals_mat.row(comp).transpose())
                                    .squaredNorm();

        auto const eidx = num_components * element_index + comp;
        _residuals->set(eidx, std::sqrt(residual / num_int_pts));
    }
}

}  // namespace NumLib
