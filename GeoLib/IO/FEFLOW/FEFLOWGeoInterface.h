/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef FEFLOWGEOINTERFACE_H_
#define FEFLOWGEOINTERFACE_H_

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

class QDomElement;
class QString;

namespace GeoLib
{
    class GEOObjects;
    class Point;
    class Polyline;
}

namespace GeoLib
{
namespace IO
{

/**
 * Read FEFLOW model files (*.fem) into OGS data structure. Currently this class supports
 * only import of mesh data and some geometry given in Supermesh section.
 */
class FEFLOWGeoInterface
{
public:
    /**
     * read a FEFLOW Model file (*.fem) in ASCII format (Version 5.4)
     *
     * This function reads mesh data in addition to geometry data given in Supermesh.
     *
     * @param filename  FEFLOW file name
     * @return a pointer to a created OGS mesh
     */
    void readFEFLOWFile(const std::string &filename, GeoLib::GEOObjects& obj);

private:
    // CLASS
    struct FEM_CLASS
    {
        unsigned problem_class = 0;
        unsigned time_mode = 0;
        unsigned orientation = 0;
        unsigned dimension = 0;
        unsigned n_layers3d = 0;
        unsigned saturation_flag = 0;
        unsigned save_fsize_rreal = 0;
        unsigned save_fsize_creal = 0;
    };

    /// read Supermesh data
    ///
    /// A super mesh is a collection of polygons, lines and points in the 2D plane
    /// and will be used for mesh generation and to define the modeling region
    void readSuperMesh(std::ifstream &feflow_file, const FEM_CLASS &fem_class, std::vector<GeoLib::Point*>** points, std::vector<GeoLib::Polyline*>** lines);

    //// read point data in Supermesh
    void readPoints(QDomElement &nodesEle, const std::string &tag, int dim, std::vector<GeoLib::Point*> &points);

};
} // IO
} // GeoLib

#endif /* FEFLOWGEOINTERFACE_H_ */
