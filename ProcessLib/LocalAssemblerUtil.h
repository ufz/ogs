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

namespace ProcessLib
{

Eigen::Map<Eigen::MatrixXd>
setupLocalMatrix(std::vector<double>& matrix_data, std::size_t const local_matrix_size);

Eigen::Map<Eigen::VectorXd>
setupLocalVector(std::vector<double>& matrix_data, std::size_t const local_matrix_size);

} // namespace ProcessLib

#endif // PROCESSLIB_LOCALASSEMBLERUTIL_H
