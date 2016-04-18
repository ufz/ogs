/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LocalAssemblerUtil.h"

namespace ProcessLib
{

Eigen::Map<Eigen::MatrixXd>
setupLocalMatrix(std::vector<double>& matrix_data, std::size_t const local_matrix_size)
{
    matrix_data.resize(local_matrix_size*local_matrix_size);
    // TODO narrowing conversion
    Eigen::Map<Eigen::MatrixXd> mat(
        matrix_data.data(), local_matrix_size, local_matrix_size);
    mat.setZero();
    return mat;
}

Eigen::Map<Eigen::VectorXd>
setupLocalVector(std::vector<double>& vector_data, std::size_t const local_vector_size)
{
    vector_data.resize(local_vector_size);
    // TODO narrowing conversion
    Eigen::Map<Eigen::VectorXd> vec(vector_data.data(), local_vector_size);
    vec.setZero();
    return vec;
}

} // namespace ProcessLib
