/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef SWMMINTERFACE_H_
#define SWMMINTERFACE_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <logog/include/logog.hpp>

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

enum class SwmmObject
{
    NODE,
    LINK,
    SUBCATCHMENT,
    SYSTEM
};

/**
 * Interface for reading files used within the
 * Storm Water Management Model (SWMM) and converting
 * those into corresponding OGS data structures.
 */
class SwmmInterface
{
public:
    SwmmInterface(std::string const& swmm_base_name);

    SwmmInterface(MeshLib::Mesh* mesh, std::string const& swmm_base_name);

    /// Reading a SWMM input file and conversion into OGS line mesh.
    MeshLib::Mesh* SwmmInputToLineMesh();

    /// Checking if file is a SWMM input file
    static int isSwmmInputFile(std::string const& inp_file_name);
    
    /// Checking if file is a SWMM output file
    int isSwmmOutputFile(std::string const& inp_file_name);

    int addResultsToMesh(MeshLib::Mesh &mesh, std::vector<double> data);

    /// Reading a SWMM input file and conversion into OGS geometry.
    static int SwmmInputToGeometry(std::string const& inp_file_name, GeoLib::GEOObjects &geo_objects);


private:
    /// Checking if file is a SWMM input file
    static int isSwmmFile(std::string const& inp_file_name, std::ifstream &in, std::string const& ext);

    /// Reading points from SWMM input file and converting them into OGS point vector.
    template <typename T>
    static int readCoordinates(std::ifstream &in, std::vector<T*> &points, std::vector<std::string> &names)
    {
        bool finished (false);
        std::size_t id (points.size());
        std::string line;
        while (!finished)
        {
            getline(in, line);
            std::size_t pos_end = 0;
            std::size_t pos_beg = line.find_first_not_of(' ', pos_end);
            pos_end = line.find_first_of(" \n", pos_beg);

            if (line.empty() || pos_beg==pos_end)
            {
                finished = true;
                return 0;
            }

            if (line.compare(pos_beg,2,";;") == 0)
                continue;

            pos_beg = line.find_first_not_of(' ', 0);
            pos_end = line.find_first_of(' ', pos_beg);
            names.push_back(line.substr(pos_beg, pos_end-pos_beg));

            T* pnt = new T(0, 0, 0, id);
            for (std::size_t i=0; i<2; ++i)
            {
                pos_beg = line.find_first_not_of(' ', pos_end);
                pos_end = line.find_first_of(' ', pos_beg);
                if (pos_beg != std::string::npos)
                    (*pnt)[i] = BaseLib::str2number<double>(line.substr(pos_beg, pos_end-pos_beg));
                else
                {
                    ERR ("SwmmInterface::readCoordinates(): error reading coordinate %d of object %s.",
                        i, names[names.size()-1].c_str());
                    return -1;
                }
            }
            points.push_back(pnt);
            id++;
        }
        return 0;
    }

    /// Read input information associated with nodes (elevation, depth, etc.)
    int readNodeData(std::ifstream &in, std::vector<MeshLib::Node*> &nodes,
        std::map<std::string, std::size_t> const& name_id_map,
        std::vector<double> &max_depth, bool read_max_depth);

    /// Read line elements
    int readLineElements(std::ifstream &in, std::vector<MeshLib::Element*> &elements,
        std::vector<MeshLib::Node*> const& nodes, std::map<std::string, std::size_t> const& name_id_map);

    /// Reading polygons from SWMM input file and converting them into OGS polyline vector.
    static int readPolygons(std::ifstream &in, std::vector<GeoLib::Polyline*> &lines,
        std::vector<std::string> &line_names, std::vector<GeoLib::Point*> &points,
        std::vector<std::string> &pnt_names);

    //variables
    std::string const _base_name;
    MeshLib::Mesh * _mesh;
};

#endif // SWMMINTERFACE_H_
