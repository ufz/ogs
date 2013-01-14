/**
 * \file
 * \author Thomas Fischer
 * \date   2010-05-03
 * \brief  Implementation of the shp to gli converter tool.
 *
 * \copyright
 * Copyright (c)  2013, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

//ShapeLib includes
#include "shapefil.h"

// STL
#include <fstream>
#include <iostream>
#include <vector>

// Base
#include "StringTools.h"

// FileIO
#include "XmlIO/XmlGmlInterface.h"
#include "XmlIO/XmlStnInterface.h"

// GEO
#include "GEOObjects.h"
#include "Point.h"
#include "ProjectData.h"
#include "Station.h"

#include "problem.h"
Problem* aproblem = NULL;

void convertPoints (DBFHandle dbf_handle,
                    std::string const& out_fname,
                    size_t x_id,
                    size_t y_id,
                    size_t z_id,
                    std::vector<size_t> const& name_component_ids,
                    std::string& points_group_name,
                    bool station)
{
	int n_records (DBFGetRecordCount (dbf_handle));
	std::cout << "writing " << n_records << " records" << std::endl;

//	out << "#POINTS" << std::endl;
//
//	for (int k(0); k<n_records; k++) {
//		double x (DBFReadDoubleAttribute( dbf_handle, k, x_id));
//		double y (DBFReadDoubleAttribute( dbf_handle, k, y_id));
//		double z (0.0);
//		if (z_id != std::numeric_limits<size_t>::max())
//			z = DBFReadDoubleAttribute( dbf_handle, k, z_id);
//		out.precision (10);
//		out.flags (std::ios::fixed);
//		out << k << " " << x << " " << y << " " << z << std::flush;
//		if (!name_component_ids.empty()) {
//			out << " $NAME ";
//			for (size_t j(0); j<name_component_ids.size(); j++) {
//				if (name_component_ids[j] != std::numeric_limits<size_t>::max()) {
//					std::string name (DBFReadStringAttribute( dbf_handle, k, name_component_ids[j]));
//					out << name.c_str() << " ";
//				}
//			}
//		}
//		out << std::endl;
//	}
//	out << "#STOP" << std::endl;

	std::vector<GeoLib::Point*>* points (new std::vector<GeoLib::Point*>);
	points->reserve (n_records);

	std::string name;
	for (int k = 0; k < n_records; k++) {
		double x(DBFReadDoubleAttribute(dbf_handle, k, x_id));
		double y(DBFReadDoubleAttribute(dbf_handle, k, y_id));
		double z(0.0);
		if (z_id != std::numeric_limits<size_t>::max()) z = DBFReadDoubleAttribute(dbf_handle, k,
						z_id);

		name = "";
		if (!name_component_ids.empty()) {
			for (size_t j(0); j < name_component_ids.size(); j++)
				if (name_component_ids[j] != std::numeric_limits<size_t>::max()) {
					name += DBFReadStringAttribute(dbf_handle, k, name_component_ids[j]);
					name += " ";
				}
		} else name = BaseLib::number2str(k);

		if (station) {
			GeoLib::Station* pnt(GeoLib::Station::createStation(name, x, y, z));
			points->push_back(pnt);
		} else {
			GeoLib::Point* pnt(new GeoLib::Point(x, y, z));
			points->push_back(pnt);
		}
	}

	GeoLib::GEOObjects* geo_objs (new GeoLib::GEOObjects());
	if (station)
		geo_objs->addStationVec(points, points_group_name);
	else
		geo_objs->addPointVec(points, points_group_name);

	std::string schema_name;
	if (station)
		schema_name = "OpenGeoSysSTN.xsd";
	else
		schema_name = "OpenGeoSysGLI.xsd";
	ProjectData* project_data (new ProjectData);
	project_data->setGEOObjects (geo_objs);
	if (station) {
		FileIO::XmlStnInterface xml (project_data, schema_name);
		xml.setNameForExport(points_group_name);
		xml.writeToFile(out_fname);
	} else {
		FileIO::XmlGmlInterface xml (project_data, schema_name);
		xml.setNameForExport(points_group_name);
		xml.writeToFile(out_fname);
	}

	delete project_data;
}

int main (int argc, char* argv[])
{
	if (argc == 1)
	{
		std::cout << "Usage: " << argv[0] << " shape_file_name" << std::endl;
		return -1;
	}

	std::string fname (argv[1]);

	/* from SHPInterface.cpp */
	int shape_type, number_of_elements;
	double padfMinBound[4], padfMaxBound[4];

	SHPHandle hSHP = SHPOpen(fname.c_str(),"rb");
	SHPGetInfo( hSHP, &number_of_elements, &shape_type, padfMinBound, padfMaxBound );

	if ((shape_type - 1) % 10 == 0)
		std::cout << "shape file contains " << number_of_elements << " points" << std::endl;
	if ( ((shape_type - 3) % 10 == 0 || (shape_type - 5) % 10 == 0))
	{
		std::cout << "shape file contains " << number_of_elements << " polylines" <<
		std::endl;
		std::cout << "this programm only handles point-input files" << std::endl;
		SHPClose(hSHP);
		return 0;
	}
	SHPClose(hSHP);
	/* end from SHPInterface */

	DBFHandle dbf_handle = DBFOpen(fname.c_str(),"rb");
	if(dbf_handle)
	{
		char* field_name (new char[256]);
		int width(0), n_decimals(0);
		size_t n_fields (DBFGetFieldCount(dbf_handle));
		std::cout << "************************************************" << std::endl;
		std::cout << "field idx | name of field | data type of field " << std::endl;
		for (size_t field_idx (0); field_idx < n_fields; field_idx++)
		{
			DBFGetFieldInfo( dbf_handle, field_idx, field_name, &width, &n_decimals);
			if (field_idx < 10)
				std::cout << "        " << field_idx << " |" << std::flush;
			else
				std::cout << "       " << field_idx << " |" << std::flush;
			std::string field_name_str (field_name);
			for (int k(0); k < (14 - (int)field_name_str.size()); k++)
				std::cout << " ";
			std::cout << field_name_str << " |" << std::flush;

			char native_field_type (DBFGetNativeFieldType  (dbf_handle, field_idx));
			switch (native_field_type)
			{
			case 'C':
				std::cout << "             string" << std::endl;
				break;
			case 'F':
				std::cout << "              float" << std::endl;
				break;
			case 'N':
				std::cout << "            numeric" << std::endl;
				break;
			default:
				std::cout << "      n_decimal " << n_decimals << std::endl;
			}
		}
		delete [] field_name;
		std::cout << "************************************************" << std::endl;

		size_t x_id, y_id, z_id;
		std::cout <<
		"please give the field idx that should be used for reading the x coordinate: " <<
		std::flush;
		std::cin >> x_id;
		std::cout <<
		"please give the field idx that should be used for reading the y coordinate: " <<
		std::flush;
		std::cin >> y_id;
		std::cout <<
		"please give the field idx that should be used for reading the z coordinate: " <<
		std::flush;
		std::cin >> z_id;
		if (z_id > n_fields)
			z_id = std::numeric_limits<size_t>::max();

		size_t n_name_components;
		std::cout << "please give the number of fields that should be added to name: " <<
		std::flush;
		std::cin >> n_name_components;
		std::vector<size_t> name_component_ids (n_name_components,
		                                        std::numeric_limits<size_t>::max());
		if (n_name_components != 0)
			for (size_t j(0); j < n_name_components; j++)
			{
				std::cout <<
				"- please give the field idx that should be used for reading the name: "
				          << std::flush;
				std::cin >> name_component_ids[j];
			}
		for (size_t j(0); j < n_name_components; j++)
			if (name_component_ids[j] > n_fields)
				name_component_ids[j] = std::numeric_limits<size_t>::max();

		size_t station (0);

		std::cout <<
		"Should I read the information as GeoLib::Station (0) or as GeoLib::Point (1)? Please give the number: "
		          << std::flush;
		std::cin >> station;

		std::string fname_base (fname);
		if (station == 0)
			fname += ".stn";
		else
			fname += ".gml";

		std::cout << "writing " << fname << " ... " << std::flush;
		convertPoints (dbf_handle,
		               fname,
		               x_id,
		               y_id,
		               z_id,
		               name_component_ids,
		               fname_base,
		               station == 0 ? true : false);
		DBFClose (dbf_handle);
		std::cout << "ok" << std::endl;
	}
	else
		std::cout << "error" << std::endl;

	return 0;
}
