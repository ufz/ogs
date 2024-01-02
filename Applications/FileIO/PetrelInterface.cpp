/**
 * \file
 * \author Thomas Fischer
 * \date   2010-02-16
 * \brief  Implementation of the PetrelInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PetrelInterface.h"

#include <fstream>

#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/StationBorehole.h"

namespace FileIO
{
PetrelInterface::PetrelInterface(std::list<std::string> const& sfc_fnames,
                                 std::list<std::string> const& well_path_fnames,
                                 std::string& unique_model_name,
                                 GeoLib::GEOObjects* geo_obj)
    : _unique_name(unique_model_name)
{
    for (auto const& surface_fname : sfc_fnames)
    {
        INFO("PetrelInterface::PetrelInterface(): open surface file.");
        std::ifstream in(surface_fname);
        if (in)
        {
            INFO("PetrelInterface::PetrelInterface(): \tdone.");
            readPetrelSurfacePoints(in);
            in.close();
        }
        else
        {
            WARN(
                "PetrelInterface::PetrelInterface(): \tCould not open file "
                "{:s}.",
                surface_fname);
        }
    }

    for (auto const& well_path_fname : well_path_fnames)
    {
        INFO("PetrelInterface::PetrelInterface(): open well path file.");
        std::ifstream in(well_path_fname);
        if (in)
        {
            INFO("PetrelInterface::PetrelInterface(): \tdone.");
            readPetrelWellTrace(in);
            in.close();
        }
        else
        {
            WARN(
                "PetrelInterface::PetrelInterface(): \tCould not open well "
                "path file {:s}.",
                well_path_fname);
        }
    }

    // move data to GEOObject
    geo_obj->addPointVec(std::move(pnt_vec), _unique_name);
    if (!well_vec.empty())
    {
        geo_obj->addStationVec(std::move(well_vec), _unique_name);
    }
}

static std::string readLine(std::istream& in)
{
    std::string line;
    std::getline(in, line);
    return line;
}

static std::list<std::string> split(std::string const& line)
{
    return BaseLib::splitString(line, ' ');
}

void PetrelInterface::readPetrelSurfacePoints(std::istream& in)
{
    std::string line = readLine(in);

    if (line.find("# Petrel Points with attributes") != std::string::npos)
    {
        // read header
        // read Version string
        readLine(in);
        // read string BEGIN HEADER
        readLine(in);

        line = readLine(in);
        while (line.find("END HEADER") == std::string::npos)
        {
            line = readLine(in);
        }

        // read points
        while (in)
        {
            auto point = std::make_unique<GeoLib::Point>();
            in >> (*point)[0] >> (*point)[1] >> (*point)[2];
            if (in)
            {
                pnt_vec.push_back(point.release());
            }
        }
    }
    else
    {
        WARN(
            "PetrelInterface::readPetrelSurface(): problem reading petrel "
            "points from line\n'{:s}'.",
            line);
    }
}

static void printListInfo(const std::list<std::string>& str_list,
                          std::string_view const prefix,
                          std::string_view const message)
{
    for (auto const& str : str_list)
    {
        INFO("PetrelInterface::{:s}(): {:s}: {:s}.", prefix, message, str);
    }
}

static double getLastNumberFromList(const std::list<std::string>& str_list)
{
    return std::stod(*(--str_list.end()));
}

void PetrelInterface::readPetrelWellTrace(std::istream& in)
{
    std::string line = readLine(in);

    if (line.find("# WELL TRACE FROM PETREL") == std::string::npos)
    {
        return;
    }

    auto printInfo = [](auto const& list, std::string_view const message)
    { printListInfo(list, "readPetrelWellTrace", message); };

    {
        // read header
        // read well name
        printInfo(split(readLine(in)), "well name");

        // read well head x coordinate
        auto str_list = split(readLine(in));
        printInfo(str_list, "well head x coord");
        double const well_head_x = getLastNumberFromList(str_list);

        // read well head y coordinate
        str_list = split(readLine(in));
        printInfo(str_list, "well head y coord");
        double const well_head_y = getLastNumberFromList(str_list);

        // read well KB
        str_list = split(readLine(in));
        printInfo(str_list, "well kb entry");
        double const well_kb = getLastNumberFromList(str_list);

        INFO("PetrelInterface::readPetrelWellTrace(): {:f}, {:f}, {:f}.",
             well_head_x,
             well_head_y,
             well_kb);
        double const depth = 0.0;
        std::string const borehole_name = "";
        int const date = 0;
        well_vec.push_back(new GeoLib::StationBorehole(
            well_head_x, well_head_y, well_kb, depth, borehole_name, date));

        // read well type
        readLine(in);
        //        std::string type(*((str_list.end())--));

        readPetrelWellTraceData(in);
    }
}

void PetrelInterface::readPetrelWellTraceData(std::istream& in)
{
    readLine(in);

    // read yet another header line
    std::string line = readLine(in);
    while (line[0] == '#')
    {
        line = readLine(in);
    }

    // read column information
    printListInfo(split(line), "readPetrelWellTraceData", "column information");

    // read points
    double md;
    double x;
    double y;
    double z;
    double tvd;
    double dx;
    double dy;
    double azim;
    double incl;
    double dls;
    line = readLine(in);
    while (in)
    {
        if (line.size() > 1 && line[0] != '#')
        {
            std::stringstream stream(line);
            stream >> md;
            stream >> x >> y >> z;
            //            pnt_vec->push_back (new GeoLib::Point (x,y,z));
            static_cast<GeoLib::StationBorehole*>(well_vec.back())
                ->addSoilLayer(x, y, z, "unknown");
            stream >> tvd >> dx >> dy >> azim >> incl >> dls;
        }
        line = readLine(in);
    }
}
}  // end namespace FileIO
