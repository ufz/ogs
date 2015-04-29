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

#include "MathLib/MathTools.h"
#include "MeshLib/Node.h"

namespace MeshLib
{

ElementCoordinatesMappingLocal::ElementCoordinatesMappingLocal(
    const Element& e,
    const CoordinateSystem &global_coords)
: _coords(global_coords)
{
    assert(e.getDimension() <= global_coords.getDimension());
    for(unsigned i = 0; i < e.getNNodes(); i++)
        _vec_nodes.push_back(new MeshLib::Node(*(e.getNode(i))));

    getRotationMatrixToGlobal(e, global_coords, _vec_nodes, _matR2global);
    rotateToLocal(_matR2global.transpose(), _vec_nodes);
}

ElementCoordinatesMappingLocal::~ElementCoordinatesMappingLocal()
{
	for (auto p : _vec_nodes) delete p;
}

void ElementCoordinatesMappingLocal::rotateToLocal(
    const RotationMatrix &matR2local,
    std::vector<MeshLib::Node*> &vec_nodes) const
{
    for (MeshLib::Node* node : vec_nodes)
        node->setCoords((matR2local* (*node)).getCoords());
}

void ElementCoordinatesMappingLocal::getRotationMatrixToGlobal(
    const Element &e,
    const CoordinateSystem &global_coords,
    const std::vector<MeshLib::Node*> &vec_nodes,
    RotationMatrix &matR) const
{
    const std::size_t global_dim = global_coords.getDimension();

    // compute R in x=R*x' where x are original coordinates and x' are local coordinates
    if (global_dim == e.getDimension()) {
        matR = RotationMatrix::Identity();
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
        GeoLib::getNewellPlane (vec_nodes, plane_normal, d);
        // compute a rotation matrix to XY
        GeoLib::computeRotationMatrixToXY(plane_normal, matR);
        // set a transposed matrix
        matR.transposeInPlace();
    }

}

} // MeshLib
