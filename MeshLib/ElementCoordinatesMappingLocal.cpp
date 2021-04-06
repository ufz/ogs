/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ElementCoordinatesMappingLocal.h"

#include <cassert>
#include <limits>

#include "GeoLib/AnalyticalGeometry.h"
#include "MathLib/MathTools.h"
#include "MathLib/Point3d.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

namespace detail
{
/// rotate points to local coordinates
void rotateToLocal(const MeshLib::RotationMatrix& matR2local,
                   std::vector<MathLib::Point3d>& points)
{
    std::transform(points.begin(), points.end(), points.begin(),
                   [&matR2local](auto const& pnt) { return matR2local * pnt; });
}

/// get a rotation matrix to the global coordinates
/// it computes R in x=R*x' where x is original coordinates and x' is local
/// coordinates
MeshLib::RotationMatrix getRotationMatrixToGlobal(
    const unsigned element_dimension,
    const unsigned global_dim,
    const std::vector<MathLib::Point3d>& points)
{
    Eigen::Matrix3d matR;
    // compute R in x=R*x' where x are original coordinates and x' are local
    // coordinates
    if (element_dimension == 1)
    {
        Eigen::Vector3d const xx =
            (Eigen::Map<Eigen::Vector3d const>(points[1].getCoords()) -
             Eigen::Map<Eigen::Vector3d const>(points[0].getCoords()))
                .normalized();
        if (global_dim == 2)
        {
            matR = GeoLib::compute2DRotationMatrixToX(xx);
        }
        else
        {
            matR = GeoLib::compute3DRotationMatrixToX(xx);
        }
        matR.transposeInPlace();
    }
    else if (global_dim == 3 && element_dimension == 2)
    {
        // get plane normal
        auto const [plane_normal, d] = GeoLib::getNewellPlane(points);
        // compute a rotation matrix to XY
        matR = GeoLib::computeRotationMatrixToXY(plane_normal);
        // set a transposed matrix
        matR.transposeInPlace();
    }
    return matR;
}
}  // namespace detail

namespace MeshLib
{
ElementCoordinatesMappingLocal::ElementCoordinatesMappingLocal(
    const Element& e, const unsigned global_dim)
    : _global_dim(global_dim), _matR2global(Eigen::Matrix3d::Identity())
{
    assert(e.getDimension() <= global_dim);
    _points.reserve(e.getNumberOfNodes());
    for (unsigned i = 0; i < e.getNumberOfNodes(); i++)
    {
        _points.emplace_back(*(e.getNode(i)));
    }

    auto const element_dim = e.getDimension();

    if (global_dim != element_dim)
    {
        _matR2global =
            detail::getRotationMatrixToGlobal(element_dim, global_dim, _points);
        detail::rotateToLocal(_matR2global.transpose(), _points);
    }
}

}  // namespace MeshLib
