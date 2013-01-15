/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHGENERATOR_H_
#define MESHGENERATOR_H_

#include "MeshLib/Mesh.h"

namespace MeshLib
{

namespace MeshGenerator
{

/**
 * generate a line mesh
 *
 * \param length
 * \param subdivision
 * \param origin_x
 * \param origin_y
 * \param origin_z
 */
Mesh* generateLineMesh(const double length, const std::size_t subdivision, const double origin_x, const double origin_y, const double origin_z);

/**
 * generate a regular quad mesh
 *
 * \param length
 * \param subdivision
 * \param origin_x
 * \param origin_y
 * \param origin_z
 */
Mesh* generateRegularQuadMesh(const double length, const std::size_t subdivision, const double origin_x, const double origin_y, const double origin_z);

}; //MeshGenerator

} //MeshLib

#endif //MESHGENERATOR_H_
