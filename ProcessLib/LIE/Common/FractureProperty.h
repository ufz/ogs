/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

#include "Utils.h"

namespace MeshLib
{
class Element;
}
namespace ProcessLib
{
template <typename T>
struct Parameter;
}
namespace ProcessLib
{
namespace LIE
{
struct FractureProperty
{
    int fracture_id = 0;
    int mat_id = 0;
    Eigen::Vector3d point_on_fracture;
    Eigen::Vector3d normal_vector;
    /// Rotation matrix from global to local coordinates
    Eigen::MatrixXd R;
    /// Initial aperture
    ProcessLib::Parameter<double> const* aperture0 = nullptr;
    ProcessLib::Parameter<double> const* specific_storage = nullptr;
    ProcessLib::Parameter<double> const* biot_coefficient = nullptr;
};

/// configure fracture property based on a fracture element assuming
/// a fracture is a straight line/flat plane
inline void setFractureProperty(unsigned dim, MeshLib::Element const& e,
                                FractureProperty& frac_prop)
{
    // 1st node is used but using other node is also possible, because
    // a fracture is not curving
    for (int j = 0; j < 3; j++)
        frac_prop.point_on_fracture[j] = e.getNode(0)->getCoords()[j];
    computeNormalVector(e, dim, frac_prop.normal_vector);
    frac_prop.R.resize(dim, dim);
    computeRotationMatrix(e, frac_prop.normal_vector, dim, frac_prop.R);
}

}  // namespace LIE
}  // namespace ProcessLib
