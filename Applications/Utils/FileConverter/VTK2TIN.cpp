/**
 * @file VTK2TIN.cpp
 * @author Norihiro Watanabe
 * @date Feb 04, 2014
 * @brief Converts a VTK mesh into a TIN file.
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// STL
#include <string>
#include <fstream>

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "LogogSimpleFormatter.h"

// FileIO
#include "FileIO/XmlIO/Boost/BoostVtuInterface.h"
#include "Legacy/MeshIO.h"

// MeshLib
#include "Mesh.h"
#include "Elements/Element.h"
#include "Node.h"

void writeTIN(const MeshLib::Mesh &mesh, const std::string &tinFile)
{
	std::ofstream os(tinFile.c_str());
	if (!os) {
		ERR("Error: cannot open a file %s", tinFile.c_str());
		return;
	}

	os.precision(std::numeric_limits<double>::digits10);
	for (size_t ie=0; ie<mesh.getNElements(); ie++)
	{
		auto e = mesh.getElement(ie);
		os << ie << " ";
		for (size_t in=0; in<3; in++)
		{
			auto nod = e->getNode(in);
			os << (*nod)[0] << " " << (*nod)[1] << " " << (*nod)[2] << " ";
		}
		os << "\n";
	}

	os.close();
}

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Converts VTK mesh into TIN file.", ' ', "0.1");
	TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
	                                     "the name of the file containing the input mesh", true,
	                                     "", "file name of input mesh");
	cmd.add(mesh_in);
	TCLAP::ValueArg<std::string> mesh_out("o", "TIN-output-file",
	                                      "the name of the file the TIN will be written to", true,
	                                      "", "file name of output TIN");
	cmd.add(mesh_out);
	cmd.parse(argc, argv);

	MeshLib::Mesh* mesh (FileIO::BoostVtuInterface::readVTUFile(mesh_in.getValue()));
	INFO("Mesh read: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());

	writeTIN(*mesh, mesh_out.getValue());

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
