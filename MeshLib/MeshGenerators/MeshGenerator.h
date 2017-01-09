/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHGENERATOR_H_
#define MESHGENERATOR_H_

#include <array>
#include <functional>
#include <string>
#include <vector>

#include "BaseLib/Subdivision.h"
#include "MathLib/Point3d.h"
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
    const MathLib::Point3d& origin = MathLib::ORIGIN);

/**
 * Generate regularly placed mesh nodes in 1D space
 * @param vec_x_coords  a vector of x coordinates
 * @param origin        coordinates of the left-most point
 * @return a vector of created mesh nodes
 */
std::vector<MeshLib::Node*> generateRegularNodes(
    const std::vector<double> &vec_x_coords,
    const MathLib::Point3d& origin = MathLib::ORIGIN);

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
    const MathLib::Point3d& origin = MathLib::ORIGIN);

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
    const MathLib::Point3d& origin = MathLib::ORIGIN);

/**
 * Generate regularly placed mesh nodes in 3D spaces
 * @param n_cells    an array of the number of cells in x,y,z directions
 * @param cell_size  an array of cell sizes in x,y,z directions
 * @param origin     coordinates of the left-most point
 * @return a vector of created mesh nodes
 */
std::vector<MeshLib::Node*> generateRegularNodes(const std::array<unsigned,3> &n_cells,
                                                 const std::array<double,3> &cell_size,
                                                 const MathLib::Point3d& origin);

/**
 * Generate an 1D Line-Element mesh. The mesh is generated in x-direction.
 *
 * \param div Subdivision operator
 * \param origin Optional mesh's origin (the left-most point) with MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateLineMesh(const BaseLib::ISubdivision &div,
                       MathLib::Point3d const& origin = MathLib::ORIGIN,
                       std::string const& mesh_name = "mesh");

/**
 * Generate an 1D Line-Element mesh. The mesh is generated in x-direction.
 *
 * \param length Mesh's length in x-direction.
 * \param subdivision Number of subdivisions.
 * \param origin Optional mesh's origin (the left-most point) with MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateLineMesh(const double length,
                       const std::size_t subdivision,
                       MathLib::Point3d const& origin = MathLib::ORIGIN,
                       std::string   const& mesh_name = "mesh");

/**
 * Generate an 1D Line-Element mesh. The mesh is generated in x-direction.
 *
 * \param n_cells Number of cells.
 * \param cell_size Length of Line elements
 * \param origin Optional mesh's origin (the left-most point) with MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateLineMesh(const unsigned n_cells,
                       const double   cell_size,
                       MathLib::Point3d const& origin = MathLib::ORIGIN,
                       std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 2D Quad-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param div_x Subdivision operator in x direction
 * \param div_y Subdivision operator in y direction
 * \param origin Optional mesh's origin (the left-most point) with MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularQuadMesh(const BaseLib::ISubdivision &div_x,
                              const BaseLib::ISubdivision &div_y,
                              MathLib::Point3d const& origin = MathLib::ORIGIN,
                              std::string const& mesh_name = "mesh");

/**
 * Generate a regular 2D Quad-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param length Mesh's dimensions in x- and y-directions.
 * \param subdivision Number of subdivisions.
 * \param origin Optional mesh's origin (the lower left corner) with
 *               MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularQuadMesh(const double length,
                              const std::size_t subdivision,
                              MathLib::Point3d const& origin = MathLib::ORIGIN,
                              std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 2D Quad-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param x_length Mesh's dimension in x-direction.
 * \param y_length Mesh's dimension in y-direction.
 * \param x_subdivision Number of subdivisions in x-direction.
 * \param y_subdivision Number of subdivisions in y-direction.
 * \param origin Optional mesh's origin (the lower left corner) with
 *               MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularQuadMesh(const double x_length,
                              const double y_length,
                              const std::size_t x_subdivision,
                              const std::size_t y_subdivision,
                              MathLib::Point3d const& origin = MathLib::ORIGIN,
                              std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 2D Quad-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param n_x_cells Number of cells in x-direction.
 * \param n_y_cells Number of cells in y-direction.
 * \param cell_size Edge length of Quad elements
 * \param origin Optional mesh's origin (the lower left corner) with MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularQuadMesh(const unsigned n_x_cells,
                              const unsigned n_y_cells,
                              const double   cell_size,
                              MathLib::Point3d const& origin = MathLib::ORIGIN,
                              std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 2D Quad-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param n_x_cells Number of cells in x-direction.
 * \param n_y_cells Number of cells in y-direction.
 * \param cell_size_x Edge length of Quad elements in x-direction
 * \param cell_size_y Edge length of Quad elements in y-direction
 * \param origin Optional mesh's origin (the lower left corner) with MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularQuadMesh(const unsigned n_x_cells,
                              const unsigned n_y_cells,
                              const double   cell_size_x,
                              const double   cell_size_y,
                              MathLib::Point3d const& origin = MathLib::ORIGIN,
                              std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 3D Hex-Element mesh.
 *
 * \param div_x Subdivision operator in x direction
 * \param div_y Subdivision operator in y direction
 * \param div_z Subdivision operator in z direction
 * \param origin Optional mesh's origin (the bottom lower left corner) with
 *               MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularHexMesh(const BaseLib::ISubdivision &div_x,
                             const BaseLib::ISubdivision &div_y,
                             const BaseLib::ISubdivision &div_z,
                             MathLib::Point3d const& origin = MathLib::ORIGIN,
                             std::string const& mesh_name = "mesh");

/**
 * Generate a regular 3D Hex-Element mesh.
 *
 * \param length      Mesh dimensions in x- and y- and z-directions.
 * \param subdivision Number of subdivisions.
 * \param origin      Optional mesh's origin (the bottom lower left corner) with MathLib::ORIGIN default.
 * \param mesh_name   Name of the new mesh.
 */
Mesh* generateRegularHexMesh(const double length,
                             const std::size_t subdivision,
                             MathLib::Point3d const& origin = MathLib::ORIGIN,
                             std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 3D Hex-Element mesh.
 *
 * \param x_length      Mesh dimension in x-direction.
 * \param y_length      Mesh dimension in y-direction.
 * \param z_length      Mesh dimension in z-direction.
 * \param x_subdivision Number of subdivisions in x-direction.
 * \param y_subdivision Number of subdivisions in y-direction.
 * \param z_subdivision Number of subdivisions in z-direction.
 * \param origin        Optional mesh's origin (the bottom lower left corner) with MathLib::ORIGIN default.
 * \param mesh_name     Name of the new mesh.
  */
Mesh* generateRegularHexMesh(const double x_length,
                             const double y_length,
                             const double z_length,
                             const std::size_t x_subdivision,
                             const std::size_t y_subdivision,
                             const std::size_t z_subdivision,
                             MathLib::Point3d const& origin = MathLib::ORIGIN,
                             std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 3D Hex-Element mesh.
 *
 * \param n_x_cells Number of cells in x-direction.
 * \param n_y_cells Number of cells in y-direction.
 * \param n_z_cells Number of cells in z-direction.
 * \param cell_size Edge length of Hex elements
 * \param origin    Optional mesh's origin (the lower left corner) with MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularHexMesh(const unsigned n_x_cells,
                             const unsigned n_y_cells,
                             const unsigned n_z_cells,
                             const double   cell_size,
                             MathLib::Point3d const& origin = MathLib::ORIGIN,
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
 * \param origin       Optional mesh's origin (the lower left corner) with MathLib::ORIGIN default.
 * \param mesh_name    Name of the new mesh.
 */
Mesh* generateRegularHexMesh(const unsigned n_x_cells,
                             const unsigned n_y_cells,
                             const unsigned n_z_cells,
                             const double   cell_size_x,
                             const double   cell_size_y,
                             const double   cell_size_z,
                             MathLib::Point3d const& origin = MathLib::ORIGIN,
                             std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 2D Triangle-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param div_x Subdivision operator in x direction
 * \param div_y Subdivision operator in y direction
 * \param origin Optional mesh's origin (the bottom lower left corner) with
 *               MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularTriMesh(const BaseLib::ISubdivision &div_x,
                             const BaseLib::ISubdivision &div_y,
                             MathLib::Point3d const& origin = MathLib::ORIGIN,
                             std::string const& mesh_name = "mesh");

/**
 * Generate a regular 2D Triangle-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param length Mesh's dimensions in x- and y-directions.
 * \param subdivision Number of subdivisions.
 * \param origin Optional mesh's origin (the lower left corner) with
 *               MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularTriMesh(const double length,
                              const std::size_t subdivision,
                              MathLib::Point3d const& origin = MathLib::ORIGIN,
                              std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 2D Triangle-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param x_length Mesh's dimension in x-direction.
 * \param y_length Mesh's dimension in y-direction.
 * \param x_subdivision Number of subdivisions in x-direction.
 * \param y_subdivision Number of subdivisions in y-direction.
 * \param origin Optional mesh's origin (the lower left corner) with
 *               MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularTriMesh(const double x_length,
                             const double y_length,
                             const std::size_t x_subdivision,
                             const std::size_t y_subdivision,
                             MathLib::Point3d const& origin = MathLib::ORIGIN,
                             std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 2D Triangle-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param n_x_cells Number of cells in x-direction.
 * \param n_y_cells Number of cells in y-direction.
 * \param cell_size Edge length of two equal sides of isosceles triangles
 * \param origin Optional mesh's origin (the lower left corner) with MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularTriMesh(const unsigned n_x_cells,
                             const unsigned n_y_cells,
                             const double   cell_size,
                             MathLib::Point3d const& origin = MathLib::ORIGIN,
                             std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 2D Triangle-Element mesh. The mesh is generated in the
 * x-y-plane.
 *
 * \param n_x_cells    Number of cells in x-direction.
 * \param n_y_cells    Number of cells in y-direction.
 * \param cell_size_x  Edge length of triangles in x-direction.
 * \param cell_size_y  Edge length of triangles in y-direction.
 * \param origin Optional mesh's origin (the lower left corner) with MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularTriMesh(const unsigned n_x_cells,
                             const unsigned n_y_cells,
                             const double   cell_size_x,
                             const double   cell_size_y,
                             MathLib::Point3d const& origin = MathLib::ORIGIN,
                             std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 3D Prism-Element mesh.
 *
 * \param x_length Mesh's dimension in x-direction.
 * \param y_length Mesh's dimension in y-direction.
 * \param z_length Mesh's dimension in z-direction.
 * \param x_subdivision Number of subdivisions in x-direction.
 * \param y_subdivision Number of subdivisions in y-direction.
 * \param z_subdivision Number of subdivisions in z-direction.
 * \param origin Optional mesh's origin (the lower left corner) with
 *               MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularPrismMesh(const double x_length,
                               const double y_length,
                               const double z_length,
                               const std::size_t x_subdivision,
                               const std::size_t y_subdivision,
                               const std::size_t z_subdivision,
                               MathLib::Point3d const& origin = MathLib::ORIGIN,
                               std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 2D Triangle-Element mesh.
 *
 * \param n_x_cells Number of cells in x-direction.
 * \param n_y_cells Number of cells in y-direction.
 * \param n_z_cells Number of cells in z-direction.
 * \param cell_size Edge length of two equal sides of isosceles triangles + height of prism
 * \param origin Optional mesh's origin (the lower left corner) with MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularPrismMesh(const unsigned n_x_cells,
                               const unsigned n_y_cells,
                               const unsigned n_z_cells,
                               const double cell_size,
                               MathLib::Point3d const& origin = MathLib::ORIGIN,
                               std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 2D Triangle-Element mesh.
 *
 * \param n_x_cells    Number of cells in x-direction.
 * \param n_y_cells    Number of cells in y-direction.
 * \param n_z_cells    Number of cells in z-direction.
 * \param cell_size_x  Edge length of triangles in x-direction.
 * \param cell_size_y  Edge length of triangles in y-direction.
 * \param cell_size_z  Edge length of triangles in z-direction.
 * \param origin Optional mesh's origin (the lower left corner) with MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularPrismMesh(const unsigned n_x_cells,
                               const unsigned n_y_cells,
                               const unsigned n_z_cells,
                               const double   cell_size_x,
                               const double   cell_size_y,
                               const double   cell_size_z,
                               MathLib::Point3d const& origin = MathLib::ORIGIN,
                               std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 3D Tet-Element mesh.
 *
 * This algorithm generates regular grid points and split each grid cell into six tetrahedrals.
 *
 * \param div_x Subdivision operator in x direction
 * \param div_y Subdivision operator in y direction
 * \param div_z Subdivision operator in z direction
 * \param origin Optional mesh's origin (the bottom lower left corner) with MathLib::ORIGIN default.
 * \param mesh_name Name of the new mesh.
 */
Mesh* generateRegularTetMesh(const BaseLib::ISubdivision &div_x,
                             const BaseLib::ISubdivision &div_y,
                             const BaseLib::ISubdivision &div_z,
                             MathLib::Point3d const& origin = MathLib::ORIGIN,
                             std::string const& mesh_name = "mesh");

/**
 * Generate a regular 3D Tet-Element mesh.
 *
 * This algorithm generates regular grid points and split each grid cell into six tetrahedrals.
 *
 * \param x_length      Mesh dimension in x-direction.
 * \param y_length      Mesh dimension in y-direction.
 * \param z_length      Mesh dimension in z-direction.
 * \param x_subdivision Number of subdivisions in x-direction.
 * \param y_subdivision Number of subdivisions in y-direction.
 * \param z_subdivision Number of subdivisions in z-direction.
 * \param origin        Optional mesh's origin (the bottom lower left corner) with MathLib::ORIGIN default.
 * \param mesh_name     Name of the new mesh.
  */
Mesh* generateRegularTetMesh(const double x_length,
                             const double y_length,
                             const double z_length,
                             const std::size_t x_subdivision,
                             const std::size_t y_subdivision,
                             const std::size_t z_subdivision,
                             MathLib::Point3d const& origin = MathLib::ORIGIN,
                             std::string   const& mesh_name = "mesh");

/**
 * Generate a regular 3D Tet-Element mesh.
 *
 * This algorithm generates regular grid points and split each grid cell into six tetrahedrals.
 *
 * \param n_x_cells    Number of cells in x-direction.
 * \param n_y_cells    Number of cells in y-direction.
 * \param n_z_cells    Number of cells in z-direction.
 * \param cell_size_x  Edge length of Tet elements in x-direction.
 * \param cell_size_y  Edge length of Tet elements in y s-direction.
 * \param cell_size_z  Edge length of Tet elements in z-direction
 * \param origin       Optional mesh's origin (the lower left corner) with MathLib::ORIGIN default.
 * \param mesh_name    Name of the new mesh.
 */
Mesh* generateRegularTetMesh(const unsigned n_x_cells,
                             const unsigned n_y_cells,
                             const unsigned n_z_cells,
                             const double   cell_size_x,
                             const double   cell_size_y,
                             const double   cell_size_z,
                             MathLib::Point3d const& origin = MathLib::ORIGIN,
                             std::string   const& mesh_name = "mesh");

/// Constructs a surface mesh approximating a surface in the 3d space given by
/// a function.
/// The surface within the xy-domain \f$[ll[0], ur[0]] \times [ll[1], ur[1]\f$
/// is described using the function \f$f(x,y)\f$.
MeshLib::Mesh*
createSurfaceMesh(std::string const& mesh_name,
    MathLib::Point3d const& ll, MathLib::Point3d const& ur,
    std::array<std::size_t, 2> const& n_steps,
    const std::function<double(double,double)>& f);

}  //MeshGenerator
} //MeshLib

#endif //MESHGENERATOR_H_
