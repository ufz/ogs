/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
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

#include <string>

#include "GeoLib/Point.h"
#include "MeshLib/Mesh.h"

namespace MeshLib
{
namespace MeshGenerator
{
/**
 * Generate an 1D Line-Element mesh. The mesh is generated in x-direction.
 *
 * \param length Mesh's length in x-direction.
 * \param subdivision Number of subdivisions.
 * \param origin Optional mesh's origin (the left-most point) with GeoLib::ORIGIN default.
 */
Mesh* generateLineMesh(const double length,
                       const std::size_t subdivision,
                       GeoLib::Point const& origin = GeoLib::ORIGIN);

/**
 * Generate an 1D Line-Element mesh. The mesh is generated in x-direction.
 *
 * \param n_cells Number of cells.
 * \param cell_size Length of Line elements
 * \param origin Optional mesh's origin (the left-most point) with GeoLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateLineMesh(const unsigned n_cells,
                       const double   cell_size,
                       GeoLib::Point const& origin = GeoLib::ORIGIN,
					   std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 2D Quad-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param length Mesh's dimensions in x- and y-directions.
 * \param subdivision Number of subdivisions.
 * \param origin Optional mesh's origin (the lower left corner) with GeoLib::ORIGIN default.
 */
Mesh* generateRegularQuadMesh(const double length,
                              const std::size_t subdivision,
                              GeoLib::Point const& origin = GeoLib::ORIGIN);

/**
 * Generate a regular 2D Quad-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param n_x_cells Number of cells in x-direction.
 * \param n_y_cells Number of cells in y-direction.
 * \param cell_size Edge length of Quad elements
 * \param origin Optional mesh's origin (the lower left corner) with GeoLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularQuadMesh(const unsigned n_x_cells,
	                          const unsigned n_y_cells,
	                          const double   cell_size,
	                          GeoLib::Point const& origin = GeoLib::ORIGIN,
	                          std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 3D Hex-Element mesh.
 *
 * \param length Mesh's dimensions in x- and y- and z-directions.
 * \param subdivision Number of subdivisions.
 * \param origin Optional mesh's origin (the bottom lower left corner) with GeoLib::ORIGIN default.
 */
Mesh* generateRegularHexMesh(const double length,
                             const std::size_t subdivision,
                             GeoLib::Point const& origin = GeoLib::ORIGIN);

/**
 * Generate a regular 3D Hex-Element mesh.
 *
 * \param n_x_cells Number of cells in x-direction.
 * \param n_y_cells Number of cells in y-direction.
 * \param n_z_cells Number of cells in z-direction.
 * \param cell_size Edge length of Hex elements
 * \param origin Optional mesh's origin (the lower left corner) with GeoLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularHexMesh(const unsigned n_x_cells,
	                         const unsigned n_y_cells,
	                         const unsigned n_z_cells,
	                         const double   cell_size,
	                         GeoLib::Point const& origin = GeoLib::ORIGIN,
	                         std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 2D Triangle-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param length Mesh's dimensions in x- and y-directions.
 * \param subdivision Number of subdivisions.
 * \param origin Optional mesh's origin (the lower left corner) with GeoLib::ORIGIN default.
 */
Mesh* generateRegularTriMesh(const double length,
                              const std::size_t subdivision,
                              GeoLib::Point const& origin = GeoLib::ORIGIN);

/**
 * Generate a regular 2D Triangle-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param n_x_cells Number of cells in x-direction.
 * \param n_y_cells Number of cells in y-direction.
 * \param cell_size Edge length of two equal sides of isosceles triangles
 * \param origin Optional mesh's origin (the lower left corner) with GeoLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularTriMesh(const unsigned n_x_cells,
	                          const unsigned n_y_cells,
	                          const double   cell_size,
	                          GeoLib::Point const& origin = GeoLib::ORIGIN,
	                          std::string   const& mesh_name = "mesh");

}  //MeshGenerator
} //MeshLib

#endif //MESHGENERATOR_H_
