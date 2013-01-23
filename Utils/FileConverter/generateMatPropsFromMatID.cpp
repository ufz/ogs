/**
 * \file
 * \author Karsten Rink
 * \date   2011-12-19
 * \brief  Implementation of the generateMatPropsFromMatID tool.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>

//#include "GridAdapter.h"
//#include "StringTools.h"

int main (int argc, char* argv[])
{
	if (argc == 1)
	{
		std::cout << "Usage: " << argv[0] << " <msh-file>" << std::endl;
		std::cout << std::endl;
		std::cout <<
		"Creates a new file for material properties and sets the material ids in the msh-file to 0."
		          << std::endl;
		return -1;
	} else {
		std::cout << "functionality has to be reimplemented in OGS-6" << std::endl;
	}

//	GridAdapter grid(argv[1]);
//	std::vector<GridAdapter::Element*> *elems = const_cast< std::vector<GridAdapter::Element*>* >(grid.getElements());
//	size_t nElems(elems->size());
//
//	// create file
//	std::string filename(argv[1]);
//	std::string name = BaseLib::extractBaseNameWithoutExtension(filename);
//
//	std::string new_matname(name + "_prop");
//	std::ofstream out_prop( new_matname.c_str(), std::ios::out );
//	if (out_prop.is_open())
//	{
//		for (size_t i=0; i<nElems; i++)
//			out_prop << i << "\t" << (*elems)[i]->material << std::endl;
//		out_prop.close();
//	}
//	else
//	{
//		std::cout << "Error: Could not create property file..." << std::endl;
//		return -1;
//	}
//
//	// set mat ids to 0 and write new msh file
//	for (size_t i=0; i<nElems; i++)
//		(*elems)[i]->material = 0;
//	const MeshLib::CFEMesh* mesh(grid.getCFEMesh());
//	std::string new_mshname(name + "_new.msh");
//	std::cout << "writing mesh to file " << new_mshname << " ... " << std::flush;
//	FileIO::OGSMeshIO mesh_io;
//	mesh_io.setMesh(mesh);
//	mesh_io.writeToFile (new_mshname);
//	std::cout << "ok" << std::endl;
//
//	std::cout << "New files \"" << new_mshname << "\" and \"" << new_matname << "\" written." << std::endl;
//	std::cout << "Conversion finished." << std::endl;
//
//	return 1;
}
