/**
 * \file
 * \author Thomas Fischer
 * \date   2010-01-14
 * \brief  Implementation of the OGSIOVer4 class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "OGSIOVer4.h"

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

#include "Applications/FileIO/Legacy/createSurface.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"
#include "GeoLib/AnalyticalGeometry.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/TINInterface.h"
#include "GeoLib/Point.h"
#include "GeoLib/PointVec.h"
#include "GeoLib/Polygon.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Station.h"
#include "GeoLib/Surface.h"
#include "GeoLib/Triangle.h"

namespace FileIO
{
namespace Legacy
{
/**************************************************************************
   GeoLib- Function: readPoints
   Aufgabe: Lesen der GLI Points und schreiben in einen Vector
   08/2005 CC Implementation
   01/2010 TF big modifications
**************************************************************************/
/** reads the points inclusive their names from input stream in
 * using the OGS-4 file format */
std::string readPoints(std::istream& in, std::vector<GeoLib::Point*>* pnt_vec,
                       bool& zero_based_indexing,
                       std::map<std::string, std::size_t>* pnt_id_name_map)
{
    std::string line;
    std::size_t cnt(0);

    std::getline(in, line);
    // geometric key words start with the hash #
    // while not found a new key word do ...
    while (line.find('#') == std::string::npos && !in.eof() && !in.fail())
    {
        // read id and point coordinates
        std::stringstream inss(line);
        std::size_t id;
        double x;
        double y;
        double z;
        inss >> id >> x >> y >> z;
        if (!inss.fail())
        {
            if (cnt == 0)
            {
                zero_based_indexing = id == 0;
            }
            pnt_vec->push_back(new GeoLib::Point(x, y, z, id));

            // read mesh density
            if (line.find("$MD") != std::string::npos)
            {
                double mesh_density;
                std::size_t pos1(line.find_first_of('M'));
                inss.str(line.substr(pos1 + 2, std::string::npos));
                inss >> mesh_density;
            }

            // read name of point
            std::size_t pos(line.find("$NAME"));
            if (pos != std::string::npos)  // OK
            {
                std::size_t end_pos((line.substr(pos + 6)).find(' '));
                if (end_pos != std::string::npos)
                {
                    (*pnt_id_name_map)[line.substr(pos + 6, end_pos)] = id;
                }
                else
                {
                    (*pnt_id_name_map)[line.substr(pos + 6)] = id;
                }
            }

            std::size_t id_pos(line.find("$ID"));
            if (id_pos != std::string::npos)
            {
                WARN(
                    "readPoints(): found tag $ID - please use tag $NAME for "
                    "reading point names in point {:d}.",
                    cnt);
            }
            cnt++;
        }
        std::getline(in, line);
    }

    return line;
}

/** reads points from a vector */
void readPolylinePointVector(const std::string& fname,
                             std::vector<GeoLib::Point*>& pnt_vec,
                             GeoLib::Polyline* ply,
                             const std::string& path,
                             std::vector<std::string>& errors)
{
    // open file
    std::ifstream in(BaseLib::joinPaths(path, fname).c_str());
    if (!in)
    {
        WARN("readPolylinePointVector(): error opening stream from {:s}",
             fname);
        errors.push_back(
            "[readPolylinePointVector] error opening stream from " + fname);
        return;
    }

    double x;
    double y;
    double z;
    while (in)
    {
        in >> x >> y >> z;
        std::size_t pnt_id(pnt_vec.size());
        pnt_vec.push_back(new GeoLib::Point(x, y, z));
        ply->addPoint(pnt_id);
    }
}

/**************************************************************************
   GeoLib-Method: Read
   Task: Read polyline data from file
   Programming:
   03/2004 CC Implementation
   09/2004 OK file path for PLY files
   07/2005 CC PLY id
   08/2005 CC parameter
   09/2005 CC itoa - convert integer to string
   01/2010 TF cleaned method from unused variables
**************************************************************************/
/** read a single Polyline from stream in into the ply_vec-vector */
std::string readPolyline(std::istream& in,
                         std::vector<GeoLib::Polyline*>& ply_vec,
                         GeoLib::PolylineVec::NameIdMap& ply_vec_names,
                         std::vector<GeoLib::Point*>& pnt_vec,
                         bool zero_based_indexing,
                         const std::vector<std::size_t>& pnt_id_map,
                         const std::string& path,
                         std::vector<std::string>& errors)
{
    std::string line;
    std::string name_of_ply;
    auto* ply(new GeoLib::Polyline(pnt_vec));
    std::size_t type = 2;  // need an initial value

    // Loop over all phases or components
    do
    {
        in >> line;
        if (line.find("$ID") != std::string::npos)
        {                // subkeyword found CC
            in >> line;  // read value
        }
        //....................................................................
        if (line.find("$NAME") != std::string::npos)  // subkeyword found
        {
            in >> line;
            name_of_ply = line;  // read value
        }
        //....................................................................
        if (line.find("$TYPE") != std::string::npos)  // subkeyword found
        {
            in >> line;  // read value
            type = static_cast<std::size_t>(strtol(line.c_str(), nullptr, 0));
        }
        //....................................................................
        if (line.find("$EPSILON") != std::string::npos)
        {                // subkeyword found
            in >> line;  // read value
        }
        //....................................................................
        if (line.find("$MAT_GROUP") != std::string::npos)
        {                // subkeyword found
            in >> line;  // read value
        }
        //....................................................................
        if (line.find("$POINTS") != std::string::npos)  // subkeyword found
        {                                               // read the point ids
            in >> line;
            if (type != 100)
            {
                while (!in.eof() && !in.fail() && !line.empty() &&
                       (line.find('#') == std::string::npos) &&
                       (line.find('$') == std::string::npos))
                {
                    auto pnt_id(BaseLib::str2number<std::size_t>(line));
                    if (!zero_based_indexing)
                    {
                        pnt_id--;  // one based indexing
                    }
                    std::size_t ply_size(ply->getNumberOfPoints());
                    if (ply_size > 0)
                    {
                        if (ply->getPointID(ply_size - 1) != pnt_id_map[pnt_id])
                        {
                            ply->addPoint(pnt_id_map[pnt_id]);
                        }
                    }
                    else
                    {
                        ply->addPoint(pnt_id_map[pnt_id]);
                    }
                    in >> line;
                }
            }
            else
            {
                WARN(
                    "readPolyline(): polyline is an arc *** reading not "
                    "implemented");
                errors.emplace_back(
                    "[readPolyline] reading polyline as an arc is not "
                    "implemented");
            }
            // empty line or the keyword or subkeyword or end of file
        }
        //....................................................................
        if (line.find("$POINT_VECTOR") !=
            std::string::npos)  // subkeyword found
        {
            in >> line;  // read file name
            line = BaseLib::joinPaths(path, line);
            readPolylinePointVector(line, pnt_vec, ply, path, errors);
        }  // subkeyword found
    } while (line.find('#') == std::string::npos && !line.empty() && in);

    if (type != 100)
    {
        ply_vec_names.insert(
            std::pair<std::string, std::size_t>(name_of_ply, ply_vec.size()));
        ply_vec.push_back(ply);
    }

    return line;
}

/**************************************************************************
   GEOLib-Function:
   Task: polyline read function
   Programming:
   03/2004 CC Implementation
   05/2004 CC Modification
   04/2005 CC Modification calculate the minimal distance between points
reference for mesh density of line element calculation 07/2005 CC read ID of
polyline 08/2005 CC parameter 01/2010 TF changed signature of function
**************************************************************************/
/** reads polylines */
std::string readPolylines(std::istream& in,
                          std::vector<GeoLib::Polyline*>& ply_vec,
                          std::map<std::string, std::size_t>& ply_vec_names,
                          std::vector<GeoLib::Point*>& pnt_vec,
                          bool zero_based_indexing,
                          const std::vector<std::size_t>& pnt_id_map,
                          const std::string& path,
                          std::vector<std::string>& errors)
{
    if (!in)
    {
        WARN("readPolylines(): input stream error.");
        return std::string("");
    }
    std::string tag("#POLYLINE");

    while (!in.eof() && !in.fail() &&
           tag.find("#POLYLINE") != std::string::npos)
    {
        tag = readPolyline(in, ply_vec, ply_vec_names, pnt_vec,
                           zero_based_indexing, pnt_id_map, path, errors);
    }

    return tag;
}

/**************************************************************************
   GeoLib-Method: readSurface
   Task: Read surface data from input stream
   Programming:
   03/2004 OK Implementation
   05/2005 OK EPSILON
   09/2005 CC file_path_base
   01/2010 TF signature modification, reimplementation
**************************************************************************/
/** read a single Surface */
std::string readSurface(std::istream& in,
                        std::vector<GeoLib::Polygon*>& polygon_vec,
                        std::vector<GeoLib::Surface*>& sfc_vec,
                        std::map<std::string, std::size_t>& sfc_names,
                        const std::vector<GeoLib::Polyline*>* const ply_vec,
                        const std::map<std::string, std::size_t>& ply_vec_names,
                        GeoLib::PointVec& pnt_vec, std::string const& path,
                        std::vector<std::string>& errors)
{
    std::string line;
    GeoLib::Surface* sfc(nullptr);

    int type(-1);
    std::string name;
    std::size_t ply_id(std::numeric_limits<std::size_t>::max());

    do
    {
        in >> line;
        if (line.find("$ID") != std::string::npos)
        {                // subkeyword found CC
            in >> line;  // read value
        }
        //....................................................................
        if (line.find("$NAME") != std::string::npos)  // subkeyword found
        {
            in >> line;  // read value
            name = line;
        }
        //....................................................................
        if (line.find("$TYPE") != std::string::npos)  // subkeyword found
        {
            in >> line;  // read value
            type = strtol(line.c_str(), nullptr, 0);
        }
        //....................................................................
        if (line.find("$EPSILON") != std::string::npos)
        {                // subkeyword found
            in >> line;  // read value
        }
        //....................................................................
        if (line.find("$TIN") != std::string::npos)  // subkeyword found
        {
            in >> line;  // read value (file name)
            std::string const file_name = BaseLib::joinPaths(path, line);
            sfc =
                GeoLib::IO::TINInterface::readTIN(file_name, pnt_vec, &errors);
        }
        //....................................................................
        if (line.find("$MAT_GROUP") != std::string::npos)
        {                // subkeyword found
            in >> line;  // read value
        }
        //....................................................................
        if (line.find("$POLYLINES") != std::string::npos)  // subkeyword found
        {  // read the name of the polyline(s)
            in >> line;
            while (!in.eof() && !in.fail() && !line.empty() &&
                   (line.find('#') == std::string::npos) &&
                   (line.find('$') == std::string::npos))
            {
                if (ply_vec == nullptr)
                {
                    OGS_FATAL("The polyline vector is not allocated.");
                }

                // we did read the name of a polyline -> search the id for
                // polyline
                auto it(ply_vec_names.find(line));
                if (it != ply_vec_names.end())
                {
                    ply_id = it->second;
                }
                else
                {
                    ply_id = ply_vec->size();
                }

                if (ply_id == ply_vec->size())
                {
                    WARN("readSurface(): polyline for surface not found!");
                    errors.emplace_back(
                        "[readSurface] polyline for surface not found!");
                }
                else
                {
                    if (type == 3)
                    {
                        WARN(
                            "readSurface(): surface type 3: flat surface with "
                            "any normal direction - reading not implemented.");
                        errors.emplace_back(
                            "[readSurface] surface type 3: flat surface with "
                            "any normal direction - reading not implemented");
                    }
                    if (type == 2)
                    {
                        WARN(
                            "readSurface(): vertical surface (type 2) - "
                            "reading not implemented");
                        errors.emplace_back(
                            "[readSurface] vertical surface (type 2) - reading "
                            "not implemented");
                    }
                }
                in >> line;
            }
            // empty line or a keyword is found
        }
    } while (line.find('#') == std::string::npos && !line.empty() && in);

    if (!name.empty())
    {
        sfc_names.insert(
            std::pair<std::string, std::size_t>(name, sfc_vec.size()));
    }

    if (sfc)
    {
        // surface create by TIN
        sfc_vec.push_back(sfc);
    }
    else
    {
        // surface created by polygon
        if (ply_id != std::numeric_limits<std::size_t>::max())
        {
            if (ply_vec == nullptr)
            {
                OGS_FATAL("The polyline vector is not allocated.");
            }

            if (ply_id != ply_vec->size())
            {
                if ((*ply_vec)[ply_id]->isClosed())
                {
                    polygon_vec.push_back(
                        new GeoLib::Polygon(*((*ply_vec)[ply_id]), true));
                }
                else
                {
                    WARN(
                        "readSurface(): cannot create surface {:s} from "
                        "polyline {:d} since polyline is not closed.",
                        name, ply_id);
                }
            }
        }
    }

    return line;
}

/**************************************************************************
   GEOLib-Method:
   Task: Surface read function
   Programming:
   03/2004 OK Implementation
   05/2004 CC Modification
   01/2010 TF changed signature of function, big modifications
**************************************************************************/
std::string readSurfaces(
    std::istream& in, std::vector<GeoLib::Surface*>& sfc_vec,
    std::map<std::string, std::size_t>& sfc_names,
    const std::vector<GeoLib::Polyline*>* const ply_vec,
    const std::map<std::string, std::size_t>& ply_vec_names,
    GeoLib::PointVec& pnt_vec, const std::string& path,
    std::vector<std::string>& errors, GeoLib::GEOObjects& geo,
    std::string const& unique_name, std::string const& gmsh_path)
{
    if (!in.good())
    {
        WARN("readSurfaces(): input stream error.");
        return std::string("");
    }
    std::string tag("#SURFACE");

    std::vector<GeoLib::Polygon*> polygon_vec;

    while (!in.eof() && !in.fail() && tag.find("#SURFACE") != std::string::npos)
    {
        std::size_t n_polygons(polygon_vec.size());
        tag = readSurface(in, polygon_vec, sfc_vec, sfc_names, ply_vec,
                          ply_vec_names, pnt_vec, path, errors);
        if (n_polygons < polygon_vec.size())
        {
            INFO("Creating a surface by triangulation of the polyline ...");
            if (FileIO::createSurface(
                    *(dynamic_cast<GeoLib::Polyline*>(polygon_vec.back())), geo,
                    unique_name, gmsh_path))
            {
                INFO("\t done");
            }
            else
            {
                WARN(
                    "\t Creating a surface by triangulation of the polyline "
                    "failed.");
            }
        }
    }
    for (auto& k : polygon_vec)
    {
        delete k;
    }

    return tag;
}

bool readGLIFileV4(const std::string& fname,
                   GeoLib::GEOObjects& geo,
                   std::string& unique_name,
                   std::vector<std::string>& errors,
                   std::string const& gmsh_path)
{
    INFO("GeoLib::readGLIFile(): open stream from file {:s}.", fname);
    std::ifstream in(fname.c_str());
    if (!in)
    {
        WARN("GeoLib::readGLIFile(): could not open file {:s}.", fname);
        errors.push_back("[readGLIFileV4] error opening stream from " + fname);
        return false;
    }
    INFO("GeoLib::readGLIFile(): \t done.");

    std::string tag;
    while (tag.find("#POINTS") == std::string::npos && !in.eof())
    {
        std::getline(in, tag);
    }

    // read names of points into vector of strings
    std::map<std::string, std::size_t> pnt_id_names_map;

    bool zero_based_idx(true);
    std::vector<GeoLib::Point*> pnt_vec;
    INFO("GeoLib::readGLIFile(): read points from stream.");
    tag = readPoints(in, &pnt_vec, zero_based_idx, &pnt_id_names_map);
    INFO("GeoLib::readGLIFile(): \t ok, {:d} points read.", pnt_vec.size());

    unique_name = BaseLib::extractBaseName(fname);
    if (!pnt_vec.empty())
    {
        geo.addPointVec(std::move(pnt_vec), unique_name,
                        std::move(pnt_id_names_map), 1e-6);
    }

    // extract path for reading external files
    const std::string path = BaseLib::extractPath(fname);

    // read names of plys into temporary string-vec
    GeoLib::PolylineVec::NameIdMap ply_names;
    std::vector<GeoLib::Polyline*> ply_vec;
    GeoLib::PointVec& point_vec(
        *const_cast<GeoLib::PointVec*>(geo.getPointVecObj(unique_name)));
    auto geo_pnt_vec(
        const_cast<std::vector<GeoLib::Point*>&>(point_vec.getVector()));
    if (tag.find("#POLYLINE") != std::string::npos && in)
    {
        INFO("GeoLib::readGLIFile(): read polylines from stream.");
        tag = readPolylines(in, ply_vec, ply_names, geo_pnt_vec, zero_based_idx,
                            geo.getPointVecObj(unique_name)->getIDMap(), path,
                            errors);
        INFO("GeoLib::readGLIFile(): \t ok, {:d} polylines read.",
             ply_vec.size());
    }
    else
    {
        INFO("GeoLib::readGLIFile(): tag #POLYLINE not found.");
    }

    if (!ply_vec.empty())
    {
        geo.addPolylineVec(std::move(ply_vec), unique_name,
                           std::move(ply_names));
    }

    // Since ply_names is a unique_ptr and is given to the GEOObject instance
    // geo it is not usable anymore. For this reason a copy is necessary.
    std::map<std::string, std::size_t> ply_names_copy;
    if (geo.getPolylineVecObj(unique_name))
    {
        ply_names_copy = std::map<std::string, std::size_t>{
            geo.getPolylineVecObj(unique_name)->getNameIDMapBegin(),
            geo.getPolylineVecObj(unique_name)->getNameIDMapEnd()};
    }

    std::vector<GeoLib::Surface*> sfc_vec;
    GeoLib::SurfaceVec::NameIdMap sfc_names;
    if (tag.find("#SURFACE") != std::string::npos && in)
    {
        INFO("GeoLib::readGLIFile(): read surfaces from stream.");

        readSurfaces(in, sfc_vec, sfc_names, geo.getPolylineVec(unique_name),
                     ply_names_copy, point_vec, path, errors, geo, unique_name,
                     gmsh_path);
        INFO("GeoLib::readGLIFile(): \tok, {:d} surfaces read.",
             sfc_vec.size());
    }
    else
    {
        INFO("GeoLib::readGLIFile(): tag #SURFACE not found.");
    }
    in.close();

    if (!sfc_vec.empty())
    {
        geo.addSurfaceVec(std::move(sfc_vec), unique_name,
                          std::move(sfc_names));
    }

    return errors.empty();
}

std::size_t writeTINSurfaces(std::ofstream& os,
                             GeoLib::SurfaceVec const* sfcs_vec,
                             std::size_t sfc_count, std::string const& path)
{
    auto const& sfcs(sfcs_vec->getVector());
    for (auto const* sfc : sfcs)
    {
        os << "#SURFACE"
           << "\n";
        std::string sfc_name;
        if (sfcs_vec->getNameOfElementByID(sfc_count, sfc_name))
        {
            os << "\t$NAME "
               << "\n"
               << "\t\t" << sfc_name << "\n";
        }
        else
        {
            os << "\t$NAME "
               << "\n"
               << "\t\t" << sfc_count << "\n";
            sfc_name = std::to_string(sfc_count);
        }
        sfc_name += ".tin";
        os << "\t$TIN"
           << "\n";
        os << "\t\t" << sfc_name << "\n";
        // create tin file
        sfc_name = BaseLib::joinPaths(path, sfc_name);
        GeoLib::IO::TINInterface::writeSurfaceAsTIN(*sfc, sfc_name);
        sfc_count++;
    }
    return sfc_count;
}

void writeGLIFileV4(const std::string& fname,
                    const std::string& geo_name,
                    const GeoLib::GEOObjects& geo)
{
    std::ofstream os(fname.c_str());
    if (auto const* const pnt_vec(geo.getPointVecObj(geo_name));
        pnt_vec != nullptr)
    {
        auto const& pnts(pnt_vec->getVector());
        const std::size_t n_pnts(pnts.size());
        INFO("GeoLib::writeGLIFileV4(): writing {:d} points to file {:s}.",
             n_pnts, fname);
        os << "#POINTS"
           << "\n";
        os.precision(std::numeric_limits<double>::digits10);
        for (std::size_t k(0); k < n_pnts; k++)
        {
            os << k << " " << *(pnts[k]);
            std::string const& pnt_name(pnt_vec->getItemNameByID(k));
            if (!pnt_name.empty())
            {
                os << " $NAME " << pnt_name;
            }
            os << "\n";
        }
    }

    if (auto const* plys_vec(geo.getPolylineVecObj(geo_name));
        plys_vec != nullptr)
    {
        auto const& plys(plys_vec->getVector());
        INFO("GeoLib::writeGLIFileV4(): {:d} polylines to file {:s}.",
             plys.size(), fname);
        for (auto ply : plys)
        {
            os << "#POLYLINE"
               << "\n";
            std::string polyline_name;
            plys_vec->getNameOfElement(ply, polyline_name);
            os << " $NAME "
               << "\n"
               << "  " << polyline_name << "\n";
            os << " $POINTS"
               << "\n";
            for (std::size_t j(0); j < ply->getNumberOfPoints(); j++)
            {
                os << "  " << ply->getPointID(j) << "\n";
            }
        }
    }

    // writing surfaces as TIN files
    if (const auto* sfcs_vec(geo.getSurfaceVecObj(geo_name));
        sfcs_vec != nullptr)
    {
        writeTINSurfaces(os, sfcs_vec, 0, BaseLib::extractPath(fname));
    }

    os << "#STOP"
       << "\n";
    os.close();
}

void writeAllDataToGLIFileV4(const std::string& fname,
                             const GeoLib::GEOObjects& geo)
{
    auto const geo_names = geo.getGeometryNames();

    // extract path for reading external files
    const std::string path = BaseLib::extractPath(fname);

    std::ofstream os(fname.c_str());

    std::size_t pnts_offset(0);
    std::vector<std::size_t> pnts_id_offset;
    pnts_id_offset.push_back(0);

    // writing all points
    os << "#POINTS"
       << "\n";
    for (auto const& geo_name : geo_names)
    {
        os.precision(std::numeric_limits<double>::digits10);
        GeoLib::PointVec const* const pnt_vec(geo.getPointVecObj(geo_name));
        if (pnt_vec)
        {
            auto const& pnts(pnt_vec->getVector());
            const std::size_t n_pnts(pnts.size());
            for (std::size_t k(0); k < n_pnts; k++)
            {
                os << pnts_offset + k << " " << *(pnts[k]);
                std::string const& pnt_name(pnt_vec->getItemNameByID(k));
                if (!pnt_name.empty())
                {
                    os << "$NAME " << pnt_name;
                }
                os << "\n";
            }
            pnts_offset += pnts.size();
            pnts_id_offset.push_back(pnts_offset);
        }
    }

    INFO("GeoLib::writeAllDataToGLIFileV4(): wrote {:d} points.", pnts_offset);

    // writing all stations
    std::vector<std::string> stn_names;
    geo.getStationVectorNames(stn_names);
    for (auto& stn_name : stn_names)
    {
        os.precision(std::numeric_limits<double>::digits10);
        const std::vector<GeoLib::Point*>* pnts(geo.getStationVec(stn_name));
        if (pnts)
        {
            for (std::size_t k(0); k < pnts->size(); k++)
            {
                os << k + pnts_offset << " " << *((*pnts)[k]) << " $NAME "
                   << static_cast<GeoLib::Station*>((*pnts)[k])->getName()
                   << "\n";
            }
            pnts_offset += pnts->size();
            pnts_id_offset.push_back(pnts_offset);
        }
    }

    std::size_t plys_cnt(0);

    // writing all polylines
    for (std::size_t j(0); j < geo_names.size(); j++)
    {
        const GeoLib::PolylineVec* plys_vec(
            geo.getPolylineVecObj(geo_names[j]));
        if (plys_vec)
        {
            auto const& plys(plys_vec->getVector());
            for (auto ply : plys)
            {
                os << "#POLYLINE"
                   << "\n";
                std::string ply_name;
                os << "  $NAME\n";
                if (plys_vec->getNameOfElementByID(plys_cnt, ply_name))
                {
                    os << "    " << ply_name << "\n";
                }
                else
                {
                    os << "    " << geo_names[j] << "-" << plys_cnt << "\n";
                }
                os << "  $POINTS"
                   << "\n";
                for (std::size_t l(0); l < ply->getNumberOfPoints(); l++)
                {
                    os << "    " << pnts_id_offset[j] + ply->getPointID(l)
                       << "\n";
                }
                plys_cnt++;
            }
        }
    }

    // writing surfaces as TIN files
    std::size_t sfcs_cnt(0);
    for (auto const& geo_name : geo_names)
    {
        const GeoLib::SurfaceVec* sfcs_vec(geo.getSurfaceVecObj(geo_name));
        if (sfcs_vec)
        {
            sfcs_cnt += writeTINSurfaces(os, sfcs_vec, sfcs_cnt, path);
        }
    }

    os << "#STOP"
       << "\n";
    os.close();
}

}  // namespace Legacy
}  // end namespace FileIO
