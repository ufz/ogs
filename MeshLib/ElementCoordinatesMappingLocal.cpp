/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
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
void rotateToLocal(
    const MeshLib::RotationMatrix &matR2local,
    std::vector<MathLib::Point3d> &points)
{
    for (auto& p : points)
        p.setCoords((matR2local*p).getCoords());
}

/// get a rotation matrix to the global coordinates
/// it computes R in x=R*x' where x is original coordinates and x' is local coordinates
void getRotationMatrixToGlobal(
    const unsigned element_dimension,
    const MeshLib::CoordinateSystem &global_coords,
    const std::vector<MathLib::Point3d> &points,
    MeshLib::RotationMatrix &matR)
{
    const std::size_t global_dim = global_coords.getDimension();

    // compute R in x=R*x' where x are original coordinates and x' are local coordinates
    if (global_dim == element_dimension) {
        matR.setIdentity();
    } else if (element_dimension == 1) {
        MathLib::Vector3 xx(points[0], points[1]);
        xx.normalize();
        if (global_dim == 2)
            GeoLib::compute2DRotationMatrixToX(xx, matR);
        else
            GeoLib::compute3DRotationMatrixToX(xx, matR);
        matR.transposeInPlace();
    } else if (global_dim == 3 && element_dimension == 2) {
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
}   // namespace detail

namespace MeshLib
{

ElementCoordinatesMappingLocal::ElementCoordinatesMappingLocal(
    const Element& e,
    const CoordinateSystem &global_coords)
: _coords(global_coords), _matR2global(3,3)
{
    assert(e.getDimension() <= global_coords.getDimension());
    _points.reserve(e.getNNodes());
    for(unsigned i = 0; i < e.getNNodes(); i++)
        _points.emplace_back(e.getNode(i)->getCoords());

    detail::getRotationMatrixToGlobal(e.getDimension(), global_coords, _points, _matR2global);
#ifdef OGS_USE_EIGEN
    detail::rotateToLocal(_matR2global.transpose(), _points);
#else
    RotationMatrix* m(_matR2global.transpose());
    detail::rotateToLocal(*m, _points);
    delete m;
#endif
}

} // MeshLib
