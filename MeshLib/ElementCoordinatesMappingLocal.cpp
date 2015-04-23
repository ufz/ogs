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
    for(size_t i = 0; i < e.getNNodes(); i++)
        _point_vec.push_back(MeshLib::Node(*(e.getNode(i))));

    getRotationMatrixToGlobal(e, global_coords, _point_vec, _matR2global);
    rotateToLocal(e, global_coords, _point_vec, _matR2global.transpose(), _point_vec);
}

void ElementCoordinatesMappingLocal::rotateToLocal(
    const Element &ele,
    const CoordinateSystem &global_coords,
    const std::vector<MeshLib::Node> &vec_pt,
    const RotationMatrix &matR2local,
    std::vector<MeshLib::Node> &local_pt) const
{
    // rotate the point coordinates
    const std::size_t global_dim = global_coords.getDimension();
    Eigen::VectorXd dx = Eigen::VectorXd::Zero(global_dim);
    Eigen::VectorXd x_new = Eigen::VectorXd::Zero(3);
    for(std::size_t i = 0; i < ele.getNNodes(); i++)
    {
        for (std::size_t j=0; j<global_dim; j++)
            dx[j] = vec_pt[i].getCoords()[j];

        x_new.head(global_dim) = matR2local * dx;
        local_pt[i] = MeshLib::Node(x_new.data());
    }
};

void ElementCoordinatesMappingLocal::getRotationMatrixToGlobal(
    const Element &e,
    const CoordinateSystem &global_coords,
    const std::vector<MeshLib::Node> &vec_pt,
    RotationMatrix &matR) const
{
    const std::size_t global_dim = global_coords.getDimension();

    // compute R in x=R*x' where x is original coordinates and x' is local coordinates
    matR = RotationMatrix::Zero(global_dim, global_dim);
    if (global_dim == e.getDimension()) {
        matR = RotationMatrix::Identity(global_dim, global_dim);
    } else if (e.getDimension() == 1) {
        MathLib::Vector3 xx(vec_pt[0], vec_pt[1]);
        xx.normalize();
        if (global_dim == 2)
            GeoLib::compute2DRotationMatrixToX(xx, matR);
        else
            GeoLib::compute3DRotationMatrixToX(xx, matR);
        matR.transposeInPlace();
    } else if (global_dim == 3 && e.getDimension() == 2) {
        std::vector<MathLib::Point3d*> pnts;
        for (auto& p : vec_pt)
            pnts.push_back(const_cast<MeshLib::Node*>(&p));

        // get plane normal
        MathLib::Vector3 plane_normal;
        double d;
        GeoLib::getNewellPlane (pnts, plane_normal, d);
        //std::cout << "pn=" << plane_normal << std::endl;
        // compute a rotation matrix to XY
        MathLib::DenseMatrix<double> matToXY(3,3,0.0);
        GeoLib::computeRotationMatrixToXY(plane_normal, matToXY);
        // set a transposed matrix
        for (size_t i=0; i<global_dim; ++i)
            for (size_t j=0; j<global_dim; ++j)
                matR(i, j) = matToXY(j,i);
    }

}

} // MeshLib
