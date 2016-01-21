/**
 * \file MoveGeometry.cpp
 * \author Karsten Rink
 * \date   2015-08-04
 * \brief  A small tool to translate geometries.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// ThirdParty
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "BaseLib/LogogSimpleFormatter.h"
#include "MathLib/Vector3.h"
#include "FileIO/XmlIO/Qt/XmlGmlInterface.h"
#include "GeoLib/GEOObjects.h"
#include <QApplication>

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);
	LOGOG_INITIALIZE();
	logog::Cout *logogCout(new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logogCout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd
		("Moves the points of a geometry by a given displacement vector.", ' ', "0.1");
	TCLAP::ValueArg<double> z_arg
		("z", "z", "displacement in z direction", false, 0.0, "z-displacement");
	cmd.add(z_arg);
	TCLAP::ValueArg<double> y_arg
		("y", "y", "displacement in y direction", false, 0.0, "y-displacement");
	cmd.add(y_arg);
	TCLAP::ValueArg<double> x_arg
		("x", "x", "displacement in x direction", false, 0.0, "x-displacement");
	cmd.add(x_arg);
	TCLAP::ValueArg<std::string> geo_output_arg
		("o","output", "output geometry file (*.gml)", true, "", "output file");
	cmd.add(geo_output_arg);
	TCLAP::ValueArg<std::string> geo_input_arg
		("i","input", "input geometry file (*.gml)", true, "", "input file");
	cmd.add(geo_input_arg);
	cmd.parse( argc, argv );

	GeoLib::GEOObjects geo_objects;
	FileIO::XmlGmlInterface xml(geo_objects);
	if (!xml.readFile(geo_input_arg.getValue()))
	{
		delete logogCout;
		delete custom_format;
		LOGOG_SHUTDOWN();
		return 1;
	}

	MathLib::Vector3 displacement(0.0, 0.0, 0.0);
	if (x_arg.isSet())
		displacement[0] = x_arg.getValue();
	if (y_arg.isSet())
		displacement[1] = y_arg.getValue();
	if (z_arg.isSet())
		displacement[2] = z_arg.getValue();

	std::vector<std::string> geo_names;
	geo_objects.getGeometryNames(geo_names);

	std::vector<GeoLib::Point*> const* point_vec = geo_objects.getPointVec(geo_names[0]);
	std::size_t const n_points = point_vec->size();
	for (std::size_t i=0; i<n_points; ++i)
		for (std::size_t c=0; c<3; ++c)
			(*(*point_vec)[i])[c] += displacement[c];

	xml.setNameForExport(geo_names[0]);
	xml.writeToFile(geo_output_arg.getValue());

	delete logogCout;
	delete custom_format;
	LOGOG_SHUTDOWN();
	return 0;
}

