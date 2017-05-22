/**
 * \file
 * \author Thomas Fischer
 * \date   2010-05-03
 * \brief  Implementation of the shp to gli converter tool.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// STL
#include <fstream>
#include <vector>

#include <tclap/CmdLine.h>

// ShapeLib
#include <shapefil.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "GeoLib/IO/XmlIO/Qt/XmlGmlInterface.h"
#include "GeoLib/IO/XmlIO/Qt/XmlStnInterface.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/Station.h"

void convertPoints (DBFHandle dbf_handle,
                    std::string const& out_fname,
                    std::size_t x_id,
                    std::size_t y_id,
                    std::size_t z_id,
                    std::vector<std::size_t> const& name_component_ids,
                    std::string& points_group_name,
                    bool station)
{
    int n_records (DBFGetRecordCount (dbf_handle));
    INFO("Reading %d records.", n_records);

    auto points = std::make_unique<std::vector<GeoLib::Point*>>();
    points->reserve (n_records);

    std::string name;
    for (int k = 0; k < n_records; k++) {
        double x(DBFReadDoubleAttribute(dbf_handle, k, x_id));
        double y(DBFReadDoubleAttribute(dbf_handle, k, y_id));
        double z(0.0);
        if (z_id != std::numeric_limits<std::size_t>::max())
            z = DBFReadDoubleAttribute(dbf_handle, k, z_id);

        name.clear();
        if (!name_component_ids.empty()) {
            for (std::size_t j(0); j < name_component_ids.size(); j++)
                if (name_component_ids[j] != std::numeric_limits<std::size_t>::max()) {
                    name += DBFReadStringAttribute(dbf_handle, k, name_component_ids[j]);
                    name += " ";
                }
        }
        else
            name = std::to_string(k);

        if (station) {
            GeoLib::Station* pnt(GeoLib::Station::createStation(name, x, y, z));
            points->push_back(pnt);
        } else {
            GeoLib::Point* pnt(new GeoLib::Point(x, y, z));
            points->push_back(pnt);
        }
    }

    GeoLib::GEOObjects geo_objs;
    if (station)
        geo_objs.addStationVec(std::move(points), points_group_name);
    else
        geo_objs.addPointVec(std::move(points), points_group_name);

    if (station) {
        GeoLib::IO::XmlStnInterface xml (geo_objs);
        xml.setNameForExport(points_group_name);
        xml.writeToFile(out_fname);
    } else {
        GeoLib::IO::XmlGmlInterface xml (geo_objs);
        xml.setNameForExport(points_group_name);
        xml.writeToFile(out_fname);
    }
}

void printFieldInformationTable(DBFHandle const& dbf_handle, std::size_t n_fields)
{
    char* field_name(new char[256]);
    int width(0), n_decimals(0);
    std::stringstream out;
    out << std::endl;
    out << "************************************************" << std::endl;
    out << "field idx | name of field | data type of field " << std::endl;
    out << "------------------------------------------------" << std::endl;
    for (std::size_t field_idx(0); field_idx < n_fields; field_idx++) {
        DBFGetFieldInfo(dbf_handle, field_idx, field_name, &width, &n_decimals);
        if (field_idx < 10)
            out << "        " << field_idx << " |";
        else
            out << "       " << field_idx << " |";
        std::string field_name_str(field_name);
        for (int k(0); k < (14 - (int) field_name_str.size()); k++)
            out << " ";
        out << field_name_str << " |";

        char native_field_type(DBFGetNativeFieldType(dbf_handle, field_idx));
        switch (native_field_type) {
        case 'C':
            out << "             string" << std::endl;
            break;
        case 'F':
            out << "              float" << std::endl;
            break;
        case 'N':
            out << "            numeric" << std::endl;
            break;
        default:
            out << "      n_decimal " << n_decimals << std::endl;
            break;
        }
    }
    delete[] field_name;
    out << "************************************************" << std::endl;
    INFO("%s", out.str().c_str());
}

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Converts points contained in shape file", ' ', "0.1");
    TCLAP::ValueArg<std::string> shapefile_arg("s",
                                               "shape-file",
                                               "the name of the shape file ",
                                               true,
                                               "",
                                               "shape file");
    cmd.add( shapefile_arg );

    cmd.parse(argc, argv);

    std::string fname (shapefile_arg.getValue());

    int shape_type, number_of_elements;
    double padfMinBound[4], padfMaxBound[4];

    SHPHandle hSHP = SHPOpen(fname.c_str(),"rb");
    if (hSHP) {
        SHPGetInfo( hSHP, &number_of_elements, &shape_type, padfMinBound, padfMaxBound );

        if ((shape_type - 1) % 10 == 0)
            INFO("Shape file contains %d points.", number_of_elements);
        if ( ((shape_type - 3) % 10 == 0 || (shape_type - 5) % 10 == 0))
        {
            ERR("Shape file contains %d polylines.", number_of_elements);
            ERR("This programm only handles only files containing points.");
            SHPClose(hSHP);
            return EXIT_SUCCESS;
        }
        SHPClose(hSHP);
    } else {
        ERR("Could not open shapefile %s.", fname.c_str());
    }

    DBFHandle dbf_handle = DBFOpen(fname.c_str(),"rb");
    if(dbf_handle)
    {
        std::size_t n_fields(DBFGetFieldCount(dbf_handle));
        printFieldInformationTable(dbf_handle, n_fields);

        std::size_t x_id, y_id, z_id;
        INFO("Please give the field idx that should be used for reading the x coordinate: ");
        std::cin >> x_id;
        INFO("Please give the field idx that should be used for reading the y coordinate: ");
        std::cin >> y_id;
        INFO("Please give the field idx that should be used for reading the z coordinate: ");
        std::cin >> z_id;

        if (z_id > n_fields)
            z_id = std::numeric_limits<std::size_t>::max();

        std::size_t n_name_components;
        INFO("Please give the number of fields that should be added to name: ");
        std::cin >> n_name_components;

        std::vector<std::size_t> name_component_ids (n_name_components,
                                                std::numeric_limits<std::size_t>::max());
        if (n_name_components != 0) {
            for (std::size_t j(0); j < n_name_components; j++)
            {
                INFO("- please give the field idx that should be used for reading the name: ");
                std::cin >> name_component_ids[j];
            }
        }
        for (std::size_t j(0); j < n_name_components; j++)
            if (name_component_ids[j] > n_fields)
                name_component_ids[j] = std::numeric_limits<std::size_t>::max();

        std::size_t station (0);

        INFO("Should I read the information as GeoLib::Station (0) or as GeoLib::Point (1)? Please give the number: ");
        std::cin >> station;

        std::string fname_base (fname);
        if (station == 0)
            fname += ".stn";
        else
            fname += ".gml";

        INFO("Writing to %s.", fname.c_str());
        convertPoints (dbf_handle,
                       fname,
                       x_id,
                       y_id,
                       z_id,
                       name_component_ids,
                       fname_base,
                       station == 0 ? true : false);
        DBFClose (dbf_handle);
        INFO("\tok.");
    } else {
        ERR("Could not open the database file.");
    }

    return EXIT_SUCCESS;
}
