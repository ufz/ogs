/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef FEFLOWGEOINTERFACE_H_
#define FEFLOWGEOINTERFACE_H_

#include <iosfwd>
#include <string>
#include <vector>

class QDomElement;

namespace GeoLib
{
class GEOObjects;
class Point;
class Polyline;
}

namespace FileIO
{
/**
 * Interface to geometric data in FEFLOW files
 */
class FEFLOWGeoInterface
{
public:
    /**
     * read a FEFLOW Model file (*.fem) in ASCII format (Version 5.4)
     *
     * This function reads geometry data given in Supermesh.
     *
     * @param filename  FEFLOW file name
     * @param geo_objects Geometric objects where imported geometry data are
     * added
     */
    void readFEFLOWFile(const std::string& filename,
                        GeoLib::GEOObjects& geo_objects);

    /// read points and polylines in Supermesh section
    ///
    /// A super mesh is a collection of polygons, lines and points in the 2D
    /// plane and will be used for mesh generation and to define the modeling
    /// region
    static void readSuperMesh(std::ifstream& feflow_file, unsigned dimension,
                              std::vector<GeoLib::Point*>*& points,
                              std::vector<GeoLib::Polyline*>*& lines);

private:
    //// read point data in Supermesh
    static void readPoints(QDomElement& nodesEle, const std::string& tag,
                           int dim, std::vector<GeoLib::Point*>& points);
};
}  // FileIO

#endif /* FEFLOWGEOINTERFACE_H_ */
