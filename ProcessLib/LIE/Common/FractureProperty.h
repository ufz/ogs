/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include <Eigen/Eigen>

#include "BranchProperty.h"
#include "JunctionProperty.h"
#include "Utils.h"

namespace MeshLib
{
class Element;
}
namespace ProcessLib
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
    ProcessLib::Parameter<double> const* aperture0 = nullptr;
    std::vector<std::unique_ptr<BranchProperty>> branches_master;
    std::vector<std::unique_ptr<BranchProperty>> branches_slave;

    virtual ~FractureProperty() = default;
};

struct FracturePropertyHM : public FractureProperty
{
    ProcessLib::Parameter<double> const* specific_storage = nullptr;
    ProcessLib::Parameter<double> const* biot_coefficient = nullptr;
};

/// configure fracture property based on a fracture element assuming
/// a fracture is a straight line/flat plane
inline void setFractureProperty(unsigned dim, MeshLib::Element const& e,
                                FractureProperty& frac_prop)
{
    // 1st node is used but using other node is also possible, because
    // a fracture is not curving
    for (int j = 0; j < 3; j++)
        frac_prop.point_on_fracture[j] = e.getCenterOfGravity().getCoords()[j];
    computeNormalVector(e, dim, frac_prop.normal_vector);
    frac_prop.R.resize(dim, dim);
    computeRotationMatrix(e, frac_prop.normal_vector, dim, frac_prop.R);
}

inline BranchProperty* createBranchProperty(MeshLib::Node const& branchNode,
                                            FractureProperty const& master_frac,
                                            FractureProperty const& slave_frac)
{
    BranchProperty* branch = new BranchProperty();
    branch->node_id = branchNode.getID();
    branch->coords = Eigen::Vector3d(branchNode.getCoords());
    branch->master_fracture_ID = master_frac.fracture_id;
    branch->slave_fracture_ID = slave_frac.fracture_id;
    // set a normal vector from the master to the slave fracture
    Eigen::Vector3d branch_vector =
        slave_frac.point_on_fracture - branch->coords;
    double sign = (branch_vector.dot(master_frac.normal_vector) < 0) ? -1 : 1;
    branch->normal_vector_branch = sign * master_frac.normal_vector;
    return branch;
}

inline JunctionProperty* createJunctionProperty(
    int junction_id,
    MeshLib::Node const& junctionNode,
    std::vector<int>
        frac_ids)
{
    JunctionProperty* junction = new JunctionProperty();
    junction->junction_id = junction_id;
    junction->node_id = junctionNode.getID();
    junction->coords = Eigen::Vector3d(junctionNode.getCoords());
    for (int j = 0; j < 2; j++)
        junction->fracture_IDs[j] = frac_ids[j];
    return junction;
}

}  // namespace LIE
}  // namespace ProcessLib
