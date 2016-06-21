/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <functional>

#include <logog/include/logog.hpp>
#include <Eigen/Core>

#include "MathLib/LinAlg/BLAS.h"
#include "NumLib/Assembler/SerialExecutor.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "NumLib/Function/Interpolation.h"
#include "LocalLinearLeastSquaresExtrapolator.h"

namespace NumLib
{

template<typename PropertyTag, typename LocalAssembler>
void
LocalLinearLeastSquaresExtrapolator<PropertyTag, LocalAssembler>::
extrapolate(LocalAssemblers const& local_assemblers, PropertyTag const property)
{
    _nodal_values.setZero();

    // counts the writes to each nodal value, i.e., the summands in order to
    // compute the average afterwards
    auto counts =
        MathLib::MatrixVectorTraits<GlobalVector>::newInstance(_nodal_values);
    counts->setZero(); // TODO BLAS?

    using Self = LocalLinearLeastSquaresExtrapolator<
        PropertyTag, LocalAssembler>;

    NumLib::SerialExecutor::executeMemberDereferenced(
        *this, &Self::extrapolateElement, local_assemblers, property, *counts);

    MathLib::BLAS::componentwiseDivide(_nodal_values, _nodal_values, *counts);
}

template<typename PropertyTag, typename LocalAssembler>
void
LocalLinearLeastSquaresExtrapolator<PropertyTag, LocalAssembler>::
calculateResiduals(LocalAssemblers const& local_assemblers,
                   PropertyTag const property)
{
    assert(static_cast<std::size_t>(_residuals.size()) == local_assemblers.size());

    using Self = LocalLinearLeastSquaresExtrapolator<
        PropertyTag, LocalAssembler>;

    NumLib::SerialExecutor::executeMemberDereferenced(
        *this, &Self::calculateResiudalElement, local_assemblers, property);
}

template<typename PropertyTag, typename LocalAssembler>
void
LocalLinearLeastSquaresExtrapolator<PropertyTag, LocalAssembler>::
extrapolateElement(std::size_t const element_index,
                   LocalAssembler const& loc_asm, PropertyTag const property,
                   GlobalVector& counts)
{
    auto const& integration_point_values = loc_asm.getIntegrationPointValues(
            property, _integration_point_values_cache);

    // number of nodes in the element
    const auto nn = loc_asm.getShapeMatrix(0).rows();
    // number of integration points in the element
    const auto ni = integration_point_values.size();

    assert(ni >= static_cast<decltype(ni)>(nn) &&
           "Least squares is not possible if there are more nodes than"
           "integration points.");

    auto& N = _local_matrix_cache; // TODO make that local?
    N.resize(nn, ni); // TODO: might reallocate very often

    for (auto int_pt=decltype(ni){0}; int_pt<ni; ++int_pt)
    {
        auto const& shp_mat = loc_asm.getShapeMatrix(int_pt);
        assert(shp_mat.rows() == nn);

        // copy shape matrix to extrapolation matrix columnwise
        N.block(0, int_pt, nn, 1) = shp_mat;
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
    Eigen::VectorXd tmp = N.transpose().householderQr().solve(integration_point_values_vec);

    // option: do the least squares computation using LDLT decomposition.
    // Eigen::VectorXd tmp =
    //         (N*N.transpose()).ldlt().solve(N*integration_point_values_vec);

    // TODO: for now always zeroth component is used
    auto const& global_indices = _local_to_global(element_index, 0).rows;

    _nodal_values.add(global_indices, tmp); // TODO does that give rise to PETSc problems?
    counts.add(global_indices, std::vector<double>(global_indices.size(), 1.0));
}

template<typename PropertyTag, typename LocalAssembler>
void
LocalLinearLeastSquaresExtrapolator<PropertyTag, LocalAssembler>::
calculateResiudalElement(std::size_t const element_index,
                         LocalAssembler const& loc_asm, PropertyTag const property)
{
    auto const& gp_vals = loc_asm.getIntegrationPointValues(
            property, _integration_point_values_cache);
    const unsigned ni = gp_vals.size();        // number of gauss points

    // TODO: for now always zeroth component is used
    const auto& global_indices = _local_to_global(element_index, 0).rows;

    // filter nodal values of the current element
    std::vector<double> nodal_vals_element;
    nodal_vals_element.resize(global_indices.size());
    for (unsigned i=0; i<global_indices.size(); ++i) {
        // TODO PETSc negative indices?
        nodal_vals_element[i] = _nodal_values[global_indices[i]];
    }

    double residual = 0.0;
    double gp_val_extrapol = 0.0;

    for (unsigned gp=0; gp<ni; ++gp)
    {
        NumLib::shapeFunctionInterpolate(
            nodal_vals_element, loc_asm.getShapeMatrix(gp), gp_val_extrapol);
        auto const& ax_m_b = gp_val_extrapol - gp_vals[gp];
        residual += ax_m_b * ax_m_b;
    }

    _residuals.set(element_index, std::sqrt(residual / ni));
}

} // namespace ProcessLib
