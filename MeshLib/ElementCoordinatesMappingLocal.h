/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>
#include <Eigen/Eigen>
#include "MathLib/Point3d.h"

namespace MeshLib
{
    class Element;
}

namespace MeshLib
{
using RotationMatrix = Eigen::Matrix<double, 3u, 3u, Eigen::RowMajor>;

/**
 * This class maps node coordinates on intrinsic coordinates of the given element.
 */
class ElementCoordinatesMappingLocal final
{
public:
    /**
     * Constructor
     * \param e          Mesh element whose node coordinates are mapped
     * \param global_dim Global dimension
     */
    ElementCoordinatesMappingLocal(const Element &e, const unsigned global_dim);

    /// return the global dimension
    unsigned getGlobalDimension() const { return global_dim_; }

    /// return mapped coordinates of the node
    MathLib::Point3d const& getMappedCoordinates(std::size_t node_id) const
    {
        return points_[node_id];
    }

    /// return a rotation matrix converting to global coordinates
    const RotationMatrix& getRotationMatrixToGlobal() const {return matR2global_;}

private:
    const unsigned global_dim_;
    std::vector<MathLib::Point3d> points_;
    RotationMatrix matR2global_;
};

}  // namespace MeshLib
