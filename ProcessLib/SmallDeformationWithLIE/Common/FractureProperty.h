/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_SMALLDEFORMATION_WITH_LIE_COMMON_FRACTUREPROPERTY_H_
#define PROCESSLIB_SMALLDEFORMATION_WITH_LIE_COMMON_FRACTUREPROPERTY_H_

#include <Eigen/Eigen>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

#include "ProcessLib/Parameter/Parameter.h"

#include "FractureProperty.h"
#include "Utils.h"


namespace ProcessLib
{
namespace SmallDeformationWithLIE
{

struct FractureProperty
{
    int mat_id = 0;
    Eigen::Vector3d point_on_fracture;
    Eigen::Vector3d normal_vector;
    /// Rotation matrix from global to local coordinates
    Eigen::MatrixXd R;
    /// Initial aperture
    ProcessLib::Parameter<double> const* aperture0 = nullptr;
};

/// configure fracture property based on a fracture element assuming
/// a fracture is a straight line/flat plane
inline void setFractureProperty(unsigned dim, MeshLib::Element const& e,
                                FractureProperty &frac_prop)
{
    // 1st node is used but using other node is also possible, because
    // a fracture is not curving
    for (int j=0; j<3; j++)
        frac_prop.point_on_fracture[j] = e.getNode(0)->getCoords()[j];
    computeNormalVector(e, frac_prop.normal_vector);
    frac_prop.R.resize(dim, dim);
    computeRotationMatrix(frac_prop.normal_vector, dim, frac_prop.R);
}

}  // namespace SmallDeformationWithLIE
}  // namespace ProcessLib

#endif // PROCESSLIB_SMALLDEFORMATION_WITH_LIE_COMMON_FRACTUREPROPERTY_H_
