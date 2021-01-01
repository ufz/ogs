/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Utils.h"

#include "MeshLib/Elements/FaceRule.h"

namespace ProcessLib
{
namespace LIE
{
// ToDo (TF) change interface
void computeNormalVector(MeshLib::Element const& e, unsigned const global_dim,
                         Eigen::Vector3d& element_normal)
{
    if (global_dim == 2)
    {
        assert(e.getGeomType() == MeshLib::MeshElemType::LINE);
        Eigen::Vector3d const v1 =
            Eigen::Map<Eigen::Vector3d const>(e.getNode(1)->getCoords()) -
            Eigen::Map<Eigen::Vector3d const>(e.getNode(0)->getCoords());
        element_normal[0] = -v1[1];
        element_normal[1] = v1[0];
        element_normal[2] = 0;  // not used in 2d but needed for normalization
        element_normal.normalize();
    }
    else if (global_dim == 3)
    {
        auto const element_normal_vector =
            MeshLib::FaceRule::getSurfaceNormal(&e).normalized();
        for (int i = 0; i < 3; ++i)
        {
            element_normal[i] = element_normal_vector[i];
        }
    }
}

void computeRotationMatrix(MeshLib::Element const& e, Eigen::Vector3d const& n,
                           unsigned const global_dim, Eigen::MatrixXd& R)
{
    if (global_dim == 2)
    {
        R.resize(2, 2);
        R << n[1], -n[0], n[0], n[1];
    }
    else if (global_dim == 3)
    {
        auto const u =
            MeshLib::FaceRule::getFirstSurfaceVector(&e).normalized();
        auto const v =
            MeshLib::FaceRule::getSecondSurfaceVector(&e).normalized();

        R.resize(3, 3);
        R << u[0], u[1], u[2], v[0], v[1], v[2], n[0], n[1], n[2];
    }
}

}  // namespace LIE
}  // namespace ProcessLib
