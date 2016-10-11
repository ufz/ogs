/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef PROCESSLIB_SMALLDEFORMATION_WITH_LIE_COMMON_UTILS_H_
#define PROCESSLIB_SMALLDEFORMATION_WITH_LIE_COMMON_UTILS_H_

#include <Eigen/Eigen>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
namespace SmallDeformationWithLIE
{

/// compute a normal vector of the given element
void computeNormalVector(MeshLib::Element const& e, Eigen::Vector3d & normal_vector);

/// compute a rotation matrix from global to local coordinates based on the given normal vector
void computeRotationMatrix(Eigen::Vector3d const& normal_vector, int dim, Eigen::MatrixXd &matR);

/// compute physical coordinates from the given shape vector, i.e. from the natural coordinates
template <typename Derived>
MathLib::Point3d computePhysicalCoordinates(MeshLib::Element const&e, Eigen::MatrixBase<Derived> const& shape)
{
    MathLib::Point3d pt;
    for (unsigned i=0; i<e.getNumberOfNodes(); i++)
    {
        MeshLib::Node const& node = *e.getNode(i);
        for (unsigned j=0; j<3; j++)
            pt[j] += shape[i]*node[j];
    }
    return pt;
}

}  // namespace SmallDeformationWithLIE
}  // namespace ProcessLib

#endif
