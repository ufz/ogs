/**
 * \file
 * \author Thomas Fischer
 * \date   2010-02-16
 * \brief  Implementation of the PetrelInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file PetrelInterface.cpp
 * @date 2010-02-16
 * @author Thomas Fischer
 */

#include "PetrelInterface.h"

#include <fstream>

#include <logog/include/logog.hpp>

#include "BaseLib/StringTools.h"
#include "GeoLib/StationBorehole.h"
#include "GeoLib/GEOObjects.h"

namespace FileIO
{
PetrelInterface::PetrelInterface(std::list<std::string> &sfc_fnames,
                                 std::list<std::string> &well_path_fnames,
                                 std::string &unique_model_name, GeoLib::GEOObjects* geo_obj) :
    _unique_name(unique_model_name), pnt_vec(new std::vector<GeoLib::Point*>),
    well_vec(new std::vector<GeoLib::Point*>), ply_vec(new std::vector<GeoLib::Polyline*>)
{
    for (std::list<std::string>::const_iterator it(sfc_fnames.begin()); it
         != sfc_fnames.end(); ++it)
    {
        INFO("PetrelInterface::PetrelInterface(): open surface file.");
        std::ifstream in((*it).c_str());
        if (in)
        {
            INFO("PetrelInterface::PetrelInterface(): \tdone.");
            readPetrelSurface(in);
            in.close();
        }
        else
            WARN("PetrelInterface::PetrelInterface(): \tCould not open file %s.",
                 it->c_str());
    }

    for (std::list<std::string>::const_iterator it(well_path_fnames.begin()); it
         != well_path_fnames.end(); ++it)
    {
        INFO("PetrelInterface::PetrelInterface(): open well path file.");
        std::ifstream in((*it).c_str());
        if (in)
        {
            INFO("PetrelInterface::PetrelInterface(): \tdone.");
            readPetrelWellTrace(in);
            in.close();
        }
        else
            WARN("PetrelInterface::PetrelInterface(): \tCould not open well path file %s.", it->c_str());
    }

    // store data in GEOObject
    geo_obj->addPointVec(std::unique_ptr<std::vector<GeoLib::Point*>>(pnt_vec),
                         _unique_name);
    if (!well_vec->empty())
        geo_obj->addStationVec(
            std::unique_ptr<std::vector<GeoLib::Point*>>(well_vec),
            _unique_name);
    if (!ply_vec->empty())
        geo_obj->addPolylineVec(
            std::unique_ptr<std::vector<GeoLib::Polyline*>>(ply_vec),
            _unique_name);
}

void PetrelInterface::readPetrelSurface(std::istream &in)
{
    char buffer[MAX_COLS_PER_ROW];
    in.getline(buffer, MAX_COLS_PER_ROW);
    std::string line(buffer);

    if (line.find("# Petrel Points with attributes") != std::string::npos)
    {
        // read header
        // read Version string
        in.getline(buffer, MAX_COLS_PER_ROW);
        line = buffer;
        // read string BEGIN HEADER
        in.getline(buffer, MAX_COLS_PER_ROW);
        line = buffer;

        in.getline(buffer, MAX_COLS_PER_ROW);
        line = buffer;
        while (line.find("END HEADER") == std::string::npos)
        {
            in.getline(buffer, MAX_COLS_PER_ROW);
            line = buffer;
        }

        // read points
        std::size_t idx(pnt_vec->size());
        while (in)
        {
            pnt_vec->push_back(new GeoLib::Point);
            in >> *((*pnt_vec)[idx]);
            if (!in)
            {
                delete (*pnt_vec)[idx];
                pnt_vec->pop_back();
            }
            else
                idx++;
        }
    } else
        WARN("PetrelInterface::readPetrelSurface(): problem reading petrel points from line\n\"%s\".",
                        line.c_str());
}

void PetrelInterface::readPetrelWellTrace(std::istream &in)
{
    char buffer[MAX_COLS_PER_ROW];
    in.getline(buffer, MAX_COLS_PER_ROW);
    std::string line(buffer);

    if (line.find("# WELL TRACE FROM PETREL") != std::string::npos)
    {
        // read header
        // read well name
        in.getline(buffer, MAX_COLS_PER_ROW);
        line = buffer;
        std::list<std::string> str_list(BaseLib::splitString(line, ' '));
        std::list<std::string>::const_iterator it(str_list.begin());
        while (it != str_list.end()) {
            INFO("PetrelInterface::readPetrelWellTrace(): well name: %s.", it->c_str());
            ++it;
        }

        // read well head x coordinate
        in.getline(buffer, MAX_COLS_PER_ROW);
        line = buffer;
        str_list = BaseLib::splitString(line, ' ');
        it = str_list.begin();
        while (it != str_list.end()) {
            INFO("PetrelInterface::readPetrelWellTrace(): well head x coord: %s.", it->c_str());
            ++it;
        }
        it = (str_list.end())--;
        --it;
        char* buf;
        double well_head_x(strtod((*it).c_str(), &buf));

        // read well head y coordinate
        in.getline(buffer, MAX_COLS_PER_ROW);
        line = buffer;
        str_list = BaseLib::splitString(line, ' ');
        it = str_list.begin();
        while (it != str_list.end()) {
            INFO("PetrelInterface::readPetrelWellTrace(): well head y coord: %s.", it->c_str());
            ++it;
        }
        it = (str_list.end())--;
        --it;
        double well_head_y(strtod((*it).c_str(), &buf));

        // read well KB
        in.getline(buffer, MAX_COLS_PER_ROW);
        line = buffer;
        str_list = BaseLib::splitString(line, ' ');
        it = str_list.begin();
        while (it != str_list.end()) {
            INFO("PetrelInterface::readPetrelWellTrace(): well kb entry: %s.", it->c_str());
            ++it;
        }
        it = (str_list.end())--;
        --it;
        double well_kb(strtod((*it).c_str(), &buf));

        INFO("PetrelInterface::readPetrelWellTrace(): %f, %f, %f.",
             well_head_x,
             well_head_y,
             well_kb);
        well_vec->push_back(new GeoLib::StationBorehole(well_head_x, well_head_y, well_kb));

        // read well type
        in.getline(buffer, MAX_COLS_PER_ROW);
//        std::string type(*((str_list.end())--));

        readPetrelWellTraceData(in);
    }
}

void PetrelInterface::readPetrelWellTraceData(std::istream &in)
{
    char buffer[MAX_COLS_PER_ROW];
    in.getline(buffer, MAX_COLS_PER_ROW);
    std::string line(buffer);

    // read yet another header line
    in.getline(buffer, MAX_COLS_PER_ROW);
    line = buffer;
    while (line[0] == '#')
    {
        in.getline(buffer, MAX_COLS_PER_ROW);
        line = buffer;
    }

    // read column information
    std::list<std::string> str_list = BaseLib::splitString(line, ' ');
    auto it = str_list.begin();
    while (it != str_list.end()) {
        INFO("PetrelInterface::readPetrelWellTraceData(): column information: %s.", it->c_str());
        ++it;
    }

    // read points
    double md, x, y, z, tvd, dx, dy, azim, incl, dls;
    in.getline(buffer, MAX_COLS_PER_ROW);
    line = buffer;
    while (in)
    {
        if (line.size() > 1 && line[0] != '#')
        {
            std::stringstream stream(line);
            stream >> md;
            stream >> x >> y >> z;
            //            pnt_vec->push_back (new GeoLib::Point (x,y,z));
            static_cast<GeoLib::StationBorehole*> ((*well_vec)[well_vec->size() - 1])->addSoilLayer(
                            x, y, z, "unknown");
            stream >> tvd >> dx >> dy >> azim >> incl >> dls;
        }
        in.getline(buffer, MAX_COLS_PER_ROW);
        line = buffer;
    }
}
} // end namespace FileIO
