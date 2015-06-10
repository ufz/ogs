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
    std::vector<MathLib::Point3d*> &vec_nodes)
{
    for (auto node : vec_nodes)
        node->setCoords((matR2local*(*node)).getCoords());
}

/// get a rotation matrix to the global coordinates
/// it computes R in x=R*x' where x is original coordinates and x' is local coordinates
void getRotationMatrixToGlobal(
    const MeshLib::Element &e,
    const MeshLib::CoordinateSystem &global_coords,
    const std::vector<MathLib::Point3d*> &vec_nodes,
    MeshLib::RotationMatrix &matR)
{
    const std::size_t global_dim = global_coords.getDimension();

    // compute R in x=R*x' where x are original coordinates and x' are local coordinates
    if (global_dim == e.getDimension()) {
        matR.setIdentity();
    } else if (e.getDimension() == 1) {
        MathLib::Vector3 xx(*vec_nodes[0], *vec_nodes[1]);
        xx.normalize();
        if (global_dim == 2)
            GeoLib::compute2DRotationMatrixToX(xx, matR);
        else
            GeoLib::compute3DRotationMatrixToX(xx, matR);
        matR.transposeInPlace();
    } else if (global_dim == 3 && e.getDimension() == 2) {
        // get plane normal
        MathLib::Vector3 plane_normal;
        double d;
        //std::tie(plane_normal, d) = GeoLib::getNewellPlane(vec_nodes);
        GeoLib::getNewellPlane(vec_nodes, plane_normal, d);

        // compute a rotation matrix to XY
        GeoLib::computeRotationMatrixToXY(plane_normal, matR);
        // set a transposed matrix
        matR.transposeInPlace();
    }

}
}   // namespace detail

namespace MeshLib
{

ElementCoordinatesMappingLocal::~ElementCoordinatesMappingLocal()
{
    for (auto p : _vec_nodes) delete p;
}

ElementCoordinatesMappingLocal::ElementCoordinatesMappingLocal(
    const Element& e,
    const CoordinateSystem &global_coords)
: _coords(global_coords), _matR2global(3,3)
{
    assert(e.getDimension() <= global_coords.getDimension());
    _vec_nodes.reserve(e.getNNodes());
    for(unsigned i = 0; i < e.getNNodes(); i++)
        _vec_nodes.push_back(new MathLib::Point3d(*static_cast<MathLib::Point3d const*>(e.getNode(i))));

    detail::getRotationMatrixToGlobal(e, global_coords, _vec_nodes, _matR2global);
#ifdef OGS_USE_EIGEN
    detail::rotateToLocal(_matR2global.transpose(), _vec_nodes);
#else
    RotationMatrix* m(_matR2global.transpose());
    detail::rotateToLocal(*m, _vec_nodes);
    delete m;
#endif
}

} // MeshLib
