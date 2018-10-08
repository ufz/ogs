/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

namespace ProcessLib
{
namespace LIE
{
struct BranchProperty
{
	int node_id;
	Eigen::Vector3d coords;
	int master_fracture_ID;
	int slave_fracture_ID;
	// unit vector normal to the master fracture in a direction to the slave
	Eigen::Vector3d normal_vector_branch;

    virtual ~BranchProperty() = default;
};


}  // namespace LIE
}  // namespace ProcessLib
