/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Utils.h"

#include <cmath>

#ifndef Q_MOC_RUN // to avoid Qt4 bug, https://bugreports.qt.io/browse/QTBUG-22829
#include <boost/math/constants/constants.hpp>
#endif

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
namespace SmallDeformationWithLIE
{

namespace
{

// Computed normal vector is oriented in the left direction of the given line element
// such that computeRotationMatrix2D() returns the indentity matrix for line elements
// parallel to a vector (1,0,0)
void computeNormalVector2D(MeshLib::Element const& e, Eigen::Vector3d & normal_vector)
{
    assert(e.getGeomType() == MeshLib::MeshElemType::LINE);

    auto v1 = (*e.getNode(1)) - (*e.getNode(0));
    normal_vector[0] = -v1[1];
    normal_vector[1] = v1[0];
    normal_vector[2] = 0.0;
    normal_vector.normalize();
}

// Compute a rotation matrix from global coordinates to local coordinates whose y axis
// should be same as the given normal vector
void computeRotationMatrix2D(Eigen::Vector3d const& n, Eigen::MatrixXd& matR)
{
    matR(0,0) = n[1];
    matR(0,1) = -n[0];
    matR(1,0) = n[0];
    matR(1,1) = n[1];
}

} // no named namespace


void computeNormalVector(MeshLib::Element const& e, Eigen::Vector3d & normal_vector)
{
    if (e.getDimension() == 1)
        computeNormalVector2D(e, normal_vector);
    else
        OGS_FATAL("2D elements are not supported in computeNormalVector()");
}

void computeRotationMatrix(Eigen::Vector3d const& normal_vector, int dim, Eigen::MatrixXd &matR)
{
    if (dim==2)
        computeRotationMatrix2D(normal_vector, matR);
    else
        OGS_FATAL("2D elements are not supported in computeNormalVector()");
}

}  // namespace SmallDeformationWithLIE
}  // namespace ProcessLib
