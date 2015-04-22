/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef ELEMENTCOORDINATESMAPPINGLOCAL_H_
#define ELEMENTCOORDINATESMAPPINGLOCAL_H_

#include <vector>

#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>
#endif

#include "MathLib/Point3d.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/CoordinateSystem.h"

namespace MeshLib
{

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EMatrix;

/**
 * This class maps node coordinates on intrinsic coordinates of the given element.
 */
class ElementCoordinatesMappingLocal
{
public:
    /**
     * Constructor
     * \param e                     Mesh element
     * \param org_coord_system      Original coordinate system
     */
    ElementCoordinatesMappingLocal(const Element* e, const CoordinateSystem &org_coord_system);
	
    /// Destructor
    virtual ~ElementCoordinatesMappingLocal() {}

    /// return mapped coordinates of the node
    const MathLib::Point3d* getMappedCoordinates(size_t node_id) const
    {
        return &_point_vec[node_id];
    }

    /// return a rotation matrix converting to orinal coordinates
    const EMatrix& getRotationMatrixToOriginal() const {return _matR2original;};

private:
    /// translate points to the origin
    void translate(std::vector<MathLib::Point3d> &point_vec, const MathLib::Point3d origin);

    /// flip points vertically or horizontally
    void flip(const CoordinateSystem &coordinate_system, std::vector<MathLib::Point3d> &vec_pt);

    /// rotate points
    void rotate(const Element &e, const CoordinateSystem &coordinate_system, std::vector<MathLib::Point3d> &vec_pt);

    /// get a rotation matrix to the original coordinates
    /// it computes R in x=R*x' where x is original coordinates and x' is local coordinates
    void getRotationMatrixToOriginal(const Element &e, const CoordinateSystem &coordinate_system, const std::vector<MathLib::Point3d> &vec_pt, EMatrix &matR2original);

private:
    std::vector<MathLib::Point3d> _point_vec;
    MathLib::Point3d _pt_translate;
    EMatrix _matR2original;
};

}

#endif // ELEMENTCOORDINATESMAPPINGLOCAL_H_
