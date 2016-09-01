/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ElementCoordinatesMappingLocal.h"

#include <limits>
#include <cassert>

#include "GeoLib/AnalyticalGeometry.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"
#include "MathLib/MathTools.h"
#include "MathLib/Point3d.h"
#include "MathLib/Vector3.h"

namespace detail
{
/// rotate points to local coordinates
void rotateToLocal(const MeshLib::RotationMatrix& matR2local,
                   std::vector<MathLib::Point3d>& points)
{
    for (auto& p : points)
        p = matR2local * p;
}

/// get a rotation matrix to the global coordinates
/// it computes R in x=R*x' where x is original coordinates and x' is local
/// coordinates
void getRotationMatrixToGlobal(const unsigned element_dimension,
                               const unsigned global_dim,
                               const std::vector<MathLib::Point3d>& points,
                               MeshLib::RotationMatrix& matR)
{
    // compute R in x=R*x' where x are original coordinates and x' are local
    // coordinates
    if (element_dimension == 1)
    {
        MathLib::Vector3 xx(points[0], points[1]);
        xx.normalize();
        if (global_dim == 2)
            GeoLib::compute2DRotationMatrixToX(xx, matR);
        else
            GeoLib::compute3DRotationMatrixToX(xx, matR);
        matR.transposeInPlace();
    }
    else if (global_dim == 3 && element_dimension == 2)
    {
        // get plane normal
        MathLib::Vector3 plane_normal;
        double d;
        std::tie(plane_normal, d) = GeoLib::getNewellPlane(points);

        // compute a rotation matrix to XY
        GeoLib::computeRotationMatrixToXY(plane_normal, matR);
        // set a transposed matrix
        matR.transposeInPlace();
    }
}
}  // namespace detail

namespace MeshLib
{
ElementCoordinatesMappingLocal::ElementCoordinatesMappingLocal(
    const Element& e, const unsigned global_dim)
    : _global_dim(global_dim), _matR2global(3, 3)
{
    assert(e.getDimension() <= global_dim);
    _points.reserve(e.getNumberOfNodes());
    for(unsigned i = 0; i < e.getNumberOfNodes(); i++)
        _points.emplace_back(*(e.getNode(i)));

    auto const element_dim = e.getDimension();

    if (global_dim == element_dim)
    {
        _matR2global.setIdentity();
        return;
    }

    detail::getRotationMatrixToGlobal(element_dim, global_dim, _points, _matR2global);
    detail::rotateToLocal(_matR2global.transpose(), _points);
}

} // MeshLib
