/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace ProcessLib
{
/// In-place transpose. A helper function for the integration point data
/// getters.
/// For time being Eigen's transposeInPlace doesn't work for non-square mapped
/// matrices.
template <int Components, typename StoreValuesFunction>
std::vector<double> transposeInPlace(
    StoreValuesFunction const& store_values_function)
{
    std::vector<double> result;
    store_values_function(result);
    MathLib::toMatrix<
        Eigen::Matrix<double, Eigen::Dynamic, Components, Eigen::RowMajor>>(
        result, result.size() / Components, Components) =
        MathLib::toMatrix<
            Eigen::Matrix<double, Components, Eigen::Dynamic, Eigen::RowMajor>>(
            result, Components, result.size() / Components)
            .transpose()
            .eval();

    return result;
}

template <int Components>
void transposeInPlace(std::vector<double>& values)
{
    MathLib::toMatrix<
        Eigen::Matrix<double, Eigen::Dynamic, Components, Eigen::RowMajor>>(
        values, values.size() / Components, Components) =
        MathLib::toMatrix<
            Eigen::Matrix<double, Components, Eigen::Dynamic, Eigen::RowMajor>>(
            values, Components, values.size() / Components)
            .transpose()
            .eval();
}

inline void transposeInPlace(std::vector<double>& values,
                             unsigned const num_components)
{
    MathLib::toMatrix<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        values, values.size() / num_components, num_components) =
        MathLib::toMatrix<Eigen::Matrix<double,
                                        Eigen::Dynamic,
                                        Eigen::Dynamic,
                                        Eigen::RowMajor>>(
            values, num_components, values.size() / num_components)
            .transpose()
            .eval();
}

}  // namespace ProcessLib
