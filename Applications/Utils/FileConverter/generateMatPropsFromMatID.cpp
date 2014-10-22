/**
 * \file
 * \author Karsten Rink
 * \date   2011-12-19
 * \brief  Implementation of the generateMatPropsFromMatID tool.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "tclap/CmdLine.h"
#include "LogogSimpleFormatter.h"
#include "FileTools.h"

// FileIO
#include "readMeshFromFile.h"
#include "FileIO/VtkIO/VtuInterface.h"

// MeshLib
#include "Mesh.h"
#include "Elements/Element.h"

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd(
	        "Creates a new file for material properties and sets the material ids in the msh-file to 0",
	        ' ',
	        "0.1");

	TCLAP::ValueArg<std::string> mesh_arg("m",
	                                          "mesh",
	                                          "the mesh to open from a file",
	                                          false,
	                                          "",
	                                          "filename for mesh input");
	cmd.add( mesh_arg );

	cmd.parse( argc, argv );

	// read mesh
	MeshLib::Mesh* mesh(FileIO::readMeshFromFile(mesh_arg.getValue()));

	std::vector<MeshLib::Element*> &elems = *(const_cast<std::vector<MeshLib::Element*>*>(&(mesh->getElements())));
	std::size_t nElems(elems.size());

	std::string name = BaseLib::extractBaseNameWithoutExtension(mesh_arg.getValue());
	// create file
	std::string new_matname(name + "_prop");
	std::ofstream out_prop( new_matname.c_str(), std::ios::out );
	if (out_prop.is_open())
	{
		for (std::size_t i=0; i<nElems; i++)
			out_prop << i << "\t" << elems[i]->getValue() << "\n";
		out_prop.close();
	}
	else
	{
		ERR("Could not create property \"%s\" file.", new_matname.c_str());
		return -1;
	}

	// set mat ids to 0 and write new msh file
	for (std::size_t i=0; i<nElems; i++)
		elems[i]->setValue(0);

	std::string new_mshname(name + "_new.vtu");
	INFO("Writing mesh to file \"%s\".", new_mshname.c_str());
	FileIO::VtuInterface mesh_io(mesh);
	mesh_io.writeToFile (new_mshname);

	INFO("New files \"%s\" and \"%s\" written.", new_mshname.c_str(), new_matname.c_str());
	std::cout << "Conversion finished." << std::endl;

	return 1;
}
