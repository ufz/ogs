 /**
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <logog/include/logog.hpp>

#include "GeoLib/Station.h"

namespace GeoLib {
    class GEOObjects;
    class Point;
    class Polyline;
}

namespace MeshLib {
    class Mesh;
    class Node;
    class Element;
}

namespace FileIO
{

/// SWMM object types
enum class SwmmObject
{
    SUBCATCHMENT = 0,
    NODE = 1,
    LINK = 2,
    SYSTEM = 3
};

/**
 * Interface for reading files used within the Storm Water Management Model (SWMM) and
 * converting the data therein into corresponding OGS data structures.
 * SWMM distinguishes four different object types (defined in the SwmmObject enum):
 *     * Subcatchments
 *     * Nodes (or junctions)
 *     * Links (or conduits)
 *     * System
 * Note that each object in a SWMM input file has a name.
 * Each of the four object types has a (different) number of parameters. There can be an arbitrary
 * number of nodes, links or subcatchments but only ever one system.
 * The interface can convert the SWMM input data into a geometry or a (line-)mesh.
 * For meshes, also output data can be added to nodes and elements.
 * The interface also provides routines for returning data as vectors as well as convenience
 * methods for outputting data into CSV files. CSV files will either contain all parameters for one
 * object at all time steps or all parameters for all objects of a given type at one timestep.
 */
class SwmmInterface final
{
public:
    /// Basic method to create the interface (containing a mesh) from a SWMM file
    /// Other non-static methods rely on this being called in the beginning.
    static std::unique_ptr<SwmmInterface> create(std::string const& file_name);

    /// If a mesh has already been created, this methods allows to add node- or link-arrays as property
    /// to that mesh. The data vectors can be created using the getArrayAtTimeStep() method.
    static bool addResultsToMesh(MeshLib::Mesh &mesh, SwmmObject const type,
        std::string const& vec_name, std::vector<double> const& data);

    /// Returns the mesh generated from SWMM file content.
    MeshLib::Mesh const& getMesh() const { return *_mesh; }

    /// Returns the name of the data array for the given object type and parameter index.
    std::string getArrayName(SwmmObject obj_type, std::size_t var_idx) const;

    /// Get all the object names for a given object type
    std::vector<std::string> getNames(SwmmObject obj_type) const;

    /// Returns the Name for the indexed object of the given type (or an empty string if an error occured).
    std::string getName(SwmmObject obj_type, std::size_t idx) const;

    /// Returns the number of objects of the given type.
    std::size_t getNumberOfObjects(SwmmObject obj_type) const;

    /// Returns the number of parameters (incl. pollutants) of the given type.
    std::size_t getNumberOfParameters(SwmmObject obj_type) const;

    /// Returns the number of time steps for the simulation results.
    std::size_t getNumberOfTimeSteps() const;

    /// Returns an array for a given variable at all nodes/links from a SWMM output file for a given time step.
    std::vector<double> getArrayAtTimeStep(SwmmObject obj_type, std::size_t time_step, std::size_t var_idx) const;

    /// Returns an array for a given variable for one specific object from a SWMM output file for all time steps.
    std::vector<double> getArrayForObject(SwmmObject obj_type, std::size_t obj_idx, std::size_t var_idx) const;

    /// Writes the node coordinates into double vectors for writing of CSV files
    bool getNodeCoordinateVectors(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z) const;

    /// Writes the inlet- and outlet IDs for links into vectors for writing of CSV files
    bool getLinkPointIds(std::vector<std::size_t> &inlets, std::vector<std::size_t> &outlets) const;

    /// Write a CSV file for all object of the given type at one time step
    bool writeCsvForTimestep(std::string const& file_name, SwmmObject obj_type, std::size_t time_step) const;

    /// Write a CSV file for one object of the given type for all time steps
    bool writeCsvForObject(std::string const& file_name, SwmmObject obj_type, std::size_t obj_idx) const;

    /// Checks if file is a SWMM input file
    static bool isSwmmInputFile(std::string const& inp_file_name);

    /// Checks if a SWMM output file exists for the current input
    bool existsSwmmOutputFile() const;

    /// Reading a SWMM input file and conversion into OGS geometry.
    static bool convertSwmmInputToGeometry(std::string const& inp_file_name,
        GeoLib::GEOObjects &geo_objects, bool add_subcatchments);

    /// Destructor
    ~SwmmInterface();

private:
    /// Constructor
    SwmmInterface(std::string const& swmm_base_name);

    /// Reading a SWMM input file and creating an OGS line mesh. This is automatically called when the object is created
    bool readSwmmInputToLineMesh();

    /// Reading points from SWMM input file and converting them into OGS point-type vector.
    /// This method is shared by geometry- and mesh conversion.
    template <typename T>
    static bool readCoordinates(std::ifstream &in, std::vector<T*> &points, std::vector<std::string> &names);

    /// Reads input information associated with nodes (elevation, depth, etc.)
    bool readNodeData(std::ifstream &in, std::vector<MeshLib::Node*> &nodes,
        std::map<std::string, std::size_t> const& name_id_map,
        std::vector<double> &max_depth, bool read_max_depth);

    /// Reads links/conduits and returns them as a vector of OGS line elements
    bool readLineElements(std::ifstream &in, std::vector<MeshLib::Element*> &elements,
        std::vector<MeshLib::Node*> const& nodes, std::map<std::string, std::size_t> const& name_id_map);

    /// Reads subcatchment information
    bool readSubcatchments(std::ifstream &in, std::map< std::string, std::size_t> const& name_id_map);

    /// Reads pollutant names and parameters
    bool readPollutants(std::ifstream &in);

    /// Returns the name of an array (or an empty string if errors occured)
    std::string getArrayName(SwmmObject obj_type, std::size_t var_idx, std::size_t n_pollutants) const;

    /// Reads the location of external rain gauge time series files.
    bool addRainGaugeTimeSeriesLocations(std::ifstream &in);

    /// Matches existing subcatchment names with subsequently read polylines marking the outlines of said subcatchments.
    bool matchSubcatchmentsWithPolygons(std::vector<GeoLib::Polyline*> const& lines, std::vector<std::string> const& names);

    /// Creates a temporary string vector containing all subcatchment names in correct order
    std::vector<std::string> getSubcatchmentNameMap() const;

    /// During geometry conversion, this adds elevation values to the existing point vector
    static bool addPointElevation(std::ifstream &in,
        std::vector<GeoLib::Point*> &points, std::map<std::string,
        std::size_t> const& name_id_map);

    /// During geometry conversion, this reads links (conduits/pumps/weirs) and converts them into polylines.
    static bool readLinksAsPolylines(std::ifstream &in,
        std::vector<GeoLib::Polyline*> &lines, std::vector<std::string> &line_names,
        std::vector<GeoLib::Point*> const& points, std::map<std::string, std::size_t> const& point_names);

    /// During geometry conversion, this reads polygones representing subcatchments
    /// from the SWMM input file and converts them into OGS polyline vector.
    static bool readPolygons(std::ifstream &in, std::vector<GeoLib::Polyline*> &lines,
        std::vector<std::string> &line_names, std::vector<GeoLib::Point*> &points,
        std::vector<std::string> &pnt_names);

    /// Checks if the given line string is empty. Empty strings mark the end of sections in a SWMM input file.
    static bool isSectionFinished(std::string const& str);

    /// Checks if the given line string is a comment line.
    static bool isCommentLine(std::string const& str);

    /// Subcatchment data structure.
    /// (SWMM stores ~20 parameters for subcatchments. Depending on relevance this struct might be extended in the future.)
    struct Subcatchment
    {
        std::string name;
        std::size_t rain_gauge;
        std::size_t outlet;
        double area;
        GeoLib::Polyline* outline;
    };

    //variables
    /// All files for a given SWMM simulation have the same base name.
    std::string const _base_name;
    /// Vector storing the names of all nodes/junctions
    std::vector<std::string> _id_nodename_map;
    /// Vector storing the names of all links/conduits
    std::vector<std::string> _id_linkname_map;
    /// Vector storing the names of all pollutants
    std::vector<std::string> _pollutant_names;
    /// Vector storing information about all subcatchments
    std::vector<Subcatchment> _subcatchments;
    /// Separate node vector containing points for defining subcatchment outlines
    std::vector<GeoLib::Point*> _subcatchment_points;
    /// Vector containing rain gauge information as well the position of external time series files
    std::vector< std::pair<GeoLib::Station, std::string> > _rain_gauges;
    /// Mesh generated from SWMM input (+ optional output data)
    std::unique_ptr<MeshLib::Mesh> _mesh;
};

} // namespace FileIO
