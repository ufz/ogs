/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHGENERATOR_H_
#define MESHGENERATOR_H_

#include <array>
#include <string>
#include <vector>

#include "BaseLib/Subdivision.h"
#include "GeoLib/Point.h"
#include "MeshLib/Mesh.h"

namespace MeshLib
{
class Node;

namespace MeshGenerator
{

/**
 * Generate regularly placed mesh nodes in 3D spaces
 * @param vec_xyz_coords  a vector of coordinates in x,y,z directions
 * @param origin          coordinates of the left-most point
 * @return a vector of created mesh nodes
 */
std::vector<MeshLib::Node*> generateRegularNodes(
    const std::vector<const std::vector<double>*> &vec_xyz_coords,
    const GeoLib::Point& origin = GeoLib::ORIGIN);

/**
 * Generate regularly placed mesh nodes in 1D space
 * @param vec_x_coords  a vector of x coordinates
 * @param origin        coordinates of the left-most point
 * @return a vector of created mesh nodes
 */
std::vector<MeshLib::Node*> generateRegularNodes(
    const std::vector<double> &vec_x_coords,
    const GeoLib::Point& origin = GeoLib::ORIGIN);

/**
 * Generate regularly placed mesh nodes in 1D space
 * @param vec_x_coords  a vector of x coordinates
 * @param vec_y_coords  a vector of y coordinates
 * @param origin        coordinates of the left-most point
 * @return a vector of created mesh nodes
 */
std::vector<MeshLib::Node*> generateRegularNodes(
    std::vector<double> &vec_x_coords,
    std::vector<double> &vec_y_coords,
    const GeoLib::Point& origin = GeoLib::ORIGIN);

/**
 * Generate regularly placed mesh nodes in 1D space
 * @param vec_x_coords  a vector of x coordinates
 * @param vec_y_coords  a vector of y coordinates
 * @param vec_z_coords  a vector of z coordinates
 * @param origin        coordinates of the left-most point
 * @return a vector of created mesh nodes
 */
std::vector<MeshLib::Node*> generateRegularNodes(
    std::vector<double> &vec_x_coords,
    std::vector<double> &vec_y_coords,
    std::vector<double> &vec_z_coords,
    const GeoLib::Point& origin = GeoLib::ORIGIN);

/**
 * Generate regularly placed mesh nodes in 3D spaces
 * @param n_cells    an array of the number of cells in x,y,z directions
 * @param cell_size  an array of cell sizes in x,y,z directions
 * @param origin     coordinates of the left-most point
 * @return a vector of created mesh nodes
 */
std::vector<MeshLib::Node*> generateRegularNodes(const std::array<unsigned,3> &n_cells,
                                                 const std::array<double,3> &cell_size,
                                                 const GeoLib::Point& origin);

/**
 * Generate an 1D Line-Element mesh. The mesh is generated in x-direction.
 *
 * \param div Subdivision operator
 * \param origin Optional mesh's origin (the left-most point) with GeoLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateLineMesh(const BaseLib::ISubdivision &div,
                       GeoLib::Point const& origin = GeoLib::ORIGIN,
                       std::string const& mesh_name = "mesh");

/**
 * Generate an 1D Line-Element mesh. The mesh is generated in x-direction.
 *
 * \param length Mesh's length in x-direction.
 * \param subdivision Number of subdivisions.
 * \param origin Optional mesh's origin (the left-most point) with GeoLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateLineMesh(const double length,
                       const std::size_t subdivision,
                       GeoLib::Point const& origin = GeoLib::ORIGIN,
                       std::string   const& mesh_name = "mesh");

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
 * \param div_x Subdivision operator in x direction
 * \param div_y Subdivision operator in y direction
 * \param origin Optional mesh's origin (the left-most point) with GeoLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularQuadMesh(const BaseLib::ISubdivision &div_x,
                              const BaseLib::ISubdivision &div_y,
                              GeoLib::Point const& origin = GeoLib::ORIGIN,
                              std::string const& mesh_name = "mesh");

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
 * Generate a regular 2D Quad-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param n_x_cells Number of cells in x-direction.
 * \param n_y_cells Number of cells in y-direction.
 * \param cell_size_x Edge length of Quad elements in x-direction
 * \param cell_size_y Edge length of Quad elements in y-direction
 * \param origin Optional mesh's origin (the lower left corner) with GeoLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularQuadMesh(const unsigned n_x_cells,
	                          const unsigned n_y_cells,
	                          const double   cell_size_x,
	                          const double   cell_size_y,
	                          GeoLib::Point const& origin = GeoLib::ORIGIN,
	                          std::string   const& mesh_name = "mesh");



/**
 * Generate a regular 3D Hex-Element mesh.
 *
 * \param div_x Subdivision operator in x direction
 * \param div_y Subdivision operator in y direction
 * \param div_z Subdivision operator in z direction
 * \param origin Optional mesh's origin (the bottom lower left corner) with GeoLib::ORIGIN default.
 */
Mesh* generateRegularHexMesh(const BaseLib::ISubdivision &div_x,
                             const BaseLib::ISubdivision &div_y,
                             const BaseLib::ISubdivision &div_z,
                             GeoLib::Point const& origin = GeoLib::ORIGIN,
                             std::string const& mesh_name = "mesh");

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
 * Generate a regular 3D Hex-Element mesh.
 *
 * \param n_x_cells    Number of cells in x-direction.
 * \param n_y_cells    Number of cells in y-direction.
 * \param n_z_cells    Number of cells in z-direction.
 * \param cell_size_x  Edge length of Hex elements in x-direction.
 * \param cell_size_y  Edge length of Hex elements in y s-direction.
 * \param cell_size_z  Edge length of Hex elements in z-direction
 * \param origin Optional mesh's origin (the lower left corner) with GeoLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularHexMesh(const unsigned n_x_cells,
	                         const unsigned n_y_cells,
	                         const unsigned n_z_cells,
	                         const double   cell_size_x,
	                         const double   cell_size_y,
	                         const double   cell_size_z,
	                         GeoLib::Point const& origin = GeoLib::ORIGIN,
	                         std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 2D Triangle-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param div_x Subdivision operator in x direction
 * \param div_y Subdivision operator in y direction
 * \param origin Optional mesh's origin (the bottom lower left corner) with GeoLib::ORIGIN default.
 */
Mesh* generateRegularTriMesh(const BaseLib::ISubdivision &div_x,
                             const BaseLib::ISubdivision &div_y,
                             GeoLib::Point const& origin = GeoLib::ORIGIN,
                             std::string const& mesh_name = "mesh");


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

/**
 * Generate a regular 2D Triangle-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param n_x_cells    Number of cells in x-direction.
 * \param n_y_cells    Number of cells in y-direction.
 * \param cell_size_x  Edge length of triangles in x-direction.
 * \param cell_size_y  Edge length of triangles in y-direction.
 * \param origin Optional mesh's origin (the lower left corner) with GeoLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularTriMesh(const unsigned n_x_cells,
	                          const unsigned n_y_cells,
	                          const double   cell_size_x,
	                          const double   cell_size_y,
	                          GeoLib::Point const& origin = GeoLib::ORIGIN,
	                          std::string   const& mesh_name = "mesh");

}  //MeshGenerator
} //MeshLib

#endif //MESHGENERATOR_H_
