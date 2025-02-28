/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Core>
#include <memory>

#include "BranchProperty.h"
#include "JunctionProperty.h"
#include "MaterialLib/FractureModels/Permeability/Permeability.h"
#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/Elements/Element.h"

namespace ParameterLib
{
template <typename T>
struct Parameter;
}

namespace ProcessLib
{
namespace LIE
{
struct FractureProperty
{
    int fracture_id = 0;
    int mat_id = 0;
    Eigen::Vector3d point_on_fracture;
    Eigen::Vector3d normal_vector;
    /// Rotation matrix from global to local coordinates
    Eigen::MatrixXd R;
    /// Initial aperture
    ParameterLib::Parameter<double> const& aperture0;
    std::vector<BranchProperty> branches_master;
    std::vector<BranchProperty> branches_slave;

    FractureProperty(int const fracture_id_, int const material_id,
                     ParameterLib::Parameter<double> const& initial_aperture)
        : fracture_id(fracture_id_),
          mat_id(material_id),
          aperture0(initial_aperture)
    {
    }

    virtual ~FractureProperty() = default;
};

/// configure fracture property based on a fracture element assuming
/// a fracture is a straight line/flat plane
inline void setFractureProperty(int const dim, MeshLib::Element const& e,
                                FractureProperty& frac_prop)
{
    auto& n = frac_prop.normal_vector;
    // 1st node is used but using other node is also possible, because
    // a fracture is not curving
    for (int j = 0; j < 3; j++)
    {
        frac_prop.point_on_fracture[j] = getCenterOfGravity(e).data()[j];
    }

    const MeshLib::ElementCoordinatesMappingLocal ele_local_coord(e, dim);

    // Global to local rotation matrix:
    Eigen::MatrixXd const global2local_rotation =
        ele_local_coord.getRotationMatrixToGlobal().transpose();
    n = global2local_rotation.row(dim - 1);

    frac_prop.R = global2local_rotation.topLeftCorner(dim, dim);

    DBUG("Normal vector of the fracture element {:d}: [{:g}, {:g}, {:g}]",
         e.getID(), n[0], n[1], n[2]);
}

inline BranchProperty createBranchProperty(MeshLib::Node const& branchNode,
                                           FractureProperty const& master_frac,
                                           FractureProperty const& slave_frac)
{
    BranchProperty branch{branchNode, master_frac.fracture_id,
                          slave_frac.fracture_id};

    // set a normal vector from the master to the slave fracture
    Eigen::Vector3d branch_vector =
        slave_frac.point_on_fracture - branch.coords;
    double sign = (branch_vector.dot(master_frac.normal_vector) < 0) ? -1 : 1;
    branch.normal_vector_branch = sign * master_frac.normal_vector;
    return branch;
}
}  // namespace LIE
}  // namespace ProcessLib
