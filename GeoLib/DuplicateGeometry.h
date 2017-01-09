/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DUPLICATEGEOMETRY_H_
#define DUPLICATEGEOMETRY_H_

#include <memory>
#include <string>
#include <vector>

namespace GeoLib
{
class GEOObjects;
class Point;
class Polyline;
class Surface;

/**
 * Creates a copy of a geometry within GEOObjects
 */
class DuplicateGeometry final
{

public:
    /**
     * Creates a copy of a geometry within GEOObjects
     * \param geo_objects The container for geometries
     * \param input_name  The geometry to be copied
     * \param output_name The name of the copy (note: this might be modified by GEOObjects)
     */
    DuplicateGeometry(GeoLib::GEOObjects &geo_objects,
        std::string const& input_name,
        std::string const& output_name);

    // Returns the (possibly modified) output name of the new geometry.
    std::string const& getFinalizedOutputName() const { return _output_name; }
    // Returns a reference to the copied point vector for modification
    std::vector<GeoLib::Point*>& getPointVectorCopy();
    // Returns a reference to the copied polyline vector for modification
    std::vector<GeoLib::Polyline*>& getPolylineVectorCopy();
    // Returns a reference to the copied surface vector for modification
    std::vector<GeoLib::Surface*>& getSurfaceVectorCopy();

private:
    // creates a deep copy of the designated geometry
    void duplicate(std::string const& input_name);

    // creates a deep copy of the polyline vector
    std::unique_ptr<std::vector<GeoLib::Polyline*>> copyPolylinesVector(
        std::vector<GeoLib::Polyline*> const& polylines) const;

    // creates a deep copy of the surface vector
    std::unique_ptr<std::vector<GeoLib::Surface*>> copySurfacesVector(
        std::vector<GeoLib::Surface*> const& surfaces) const;

    std::string _output_name;
    GeoLib::GEOObjects& _geo_objects;

}; // class

} // namespace

#endif /* DUPLICATEGEOMETRY_H_ */
