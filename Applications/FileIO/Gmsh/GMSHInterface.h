/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GMSHINTERFACE_H_
#define GMSHINTERFACE_H_

#include <iosfwd>

#include "BaseLib/IO/Writer.h"

#include "MathLib/LinAlg/Dense/DenseMatrix.h"

namespace GeoLib
{
class GEOObjects;
class Polygon;
}

namespace FileIO
{

namespace GMSH
{

class GMSHPoint;
class GMSHPolygonTree;
class GMSHMeshDensityStrategy;

enum class MeshDensityAlgorithm {
    FixedMeshDensity, //!< set the parameter with a fixed value
    AdaptiveMeshDensity //!< computing the mesh density employing a QuadTree
};

/**
 * \brief Reads and writes GMSH-files to and from OGS data structures.
 */
class GMSHInterface final : public BaseLib::IO::Writer
{
public:
    /**
     * @param geo_objs reference to instance of class GEOObject that maintains
     * the geometries.
     *     The instance is used for preparation geometries for writing them to
     * the gmsh file format.
     * @param include_stations_as_constraints switch to enable writing stations
     * as constraints
     * @param mesh_density_algorithm one of the mesh density algorithms (\@see
     * enum MeshDensityAlgorithm)
     * @param param1 parameter that can be used for the mesh density algorithm
     * @param param2 parameter that can be used for the mesh density algorithm
     * @param param3 parameter that can be used for the mesh density algorithm
     * @param selected_geometries vector of names of geometries, that should be
     * employed for mesh generation.
     * @param rotate if the value of the parameter is true then the input points
     * will be rotated on the \f$x\f$-\f$y\f$-plane, else the input points will
     * be (orthogonal) projected to the \f$x\f$-\f$y\f$-plane.
     * @param keep_preprocessed_geometry keep the pre-processed geometry, useful
     * for debugging the mesh creation
     */
    GMSHInterface(GeoLib::GEOObjects& geo_objs,
                  bool include_stations_as_constraints,
                  GMSH::MeshDensityAlgorithm mesh_density_algorithm,
                  double param1, double param2, std::size_t param3,
                  std::vector<std::string>& selected_geometries,
                  bool rotate = false, bool keep_preprocessed_geometry = false);

    GMSHInterface(GMSHInterface const&) = delete;
    GMSHInterface(GMSHInterface &&) = delete;
    GMSHInterface& operator=(GMSHInterface const&) = delete;
    GMSHInterface& operator=(GMSHInterface &&) = delete;

    ~GMSHInterface();

protected:
    bool write();

private:
    /**
     * 1. get and merge data from _geo_objs
     * 2. compute topological hierarchy
     * @param out
     * @todo activate error codes and hand them on to the Writer class,
     * i.e. 0 = okay, 1 = geo_objects is empty, 2 = error while merging,
     * 3 = error writing file
     */
    int writeGMSHInputFile(std::ostream & out);

    void writePoints(std::ostream& out) const;

    std::size_t _n_lines;
    std::size_t _n_plane_sfc;

    GeoLib::GEOObjects & _geo_objs;
    std::vector<std::string>& _selected_geometries;
    std::string _gmsh_geo_name;
    std::list<GMSH::GMSHPolygonTree*> _polygon_tree_list;

    std::vector<GMSH::GMSHPoint*> _gmsh_pnts;

    GMSH::GMSHMeshDensityStrategy *_mesh_density_strategy;
    /// Holds the inverse rotation matrix. The matrix is used in writePoints() to
    /// revert the rotation done in writeGMSHInputFile().
    MathLib::DenseMatrix<double> _inverse_rot_mat =
        MathLib::DenseMatrix<double>(3, 3, 0);
    /// Signals if the input points should be rotated or projected to the
    /// \f$x\f$-\f$y\f$-plane
    bool _rotate = false;
    bool _keep_preprocessed_geometry = true;
};
} // end namespace GMSH
} // end namespace FileIO

#endif /* GMSHINTERFACE_H_ */
