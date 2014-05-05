/**
 * @file   removeMeshElements.cpp
 * @author Norihiro Watanabe
 * @date   2013/10/15
 * @brief  Remove mesh elements
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "LogogSimpleFormatter.h"

// FileIO
#include "Legacy/MeshIO.h"
#include "readMeshFromFile.h"

// MeshLib
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "MeshEnums.h"
#include "MeshEditing/ElementExtraction.h"

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Remove mesh elements.", ' ', "0.1");
	TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
	                                     "the name of the file containing the input mesh", true,
	                                     "", "file name of input mesh");
	cmd.add(mesh_in);
	TCLAP::ValueArg<std::string> mesh_out("o", "mesh-output-file",
	                                      "the name of the file the mesh will be written to", true,
	                                      "", "file name of output mesh");
	cmd.add(mesh_out);
	TCLAP::SwitchArg zveArg("z", "zero-volume", "remove zero volume elements", false);
	cmd.add(zveArg);
	TCLAP::MultiArg<std::string> eleTypeArg("t", "element-type",
	                                      "element type to be removed", false, "element type");
	cmd.add(eleTypeArg);
	TCLAP::MultiArg<unsigned> matIDArg("m", "material-id",
	                                      "material id", false, "material id");
	cmd.add(matIDArg);
	cmd.parse(argc, argv);

	MeshLib::Mesh const*const mesh (FileIO::readMeshFromFile(mesh_in.getValue()));
	INFO("Mesh read: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());
	MeshLib::ElementExtraction ex(*mesh);

	// search elements IDs to be removed
	if (zveArg.isSet()) {
		const std::size_t n_removed_elements = ex.searchByZeroContent();
		INFO("%d zero volume elements found.", n_removed_elements);
	}
	if (eleTypeArg.isSet()) {
		const std::vector<std::string> eleTypeNames = eleTypeArg.getValue();
		for (auto typeName : eleTypeNames) {
			const MeshElemType type = String2MeshElemType(typeName);
			if (type == MeshElemType::INVALID) continue;
			const std::size_t n_removed_elements = ex.searchByElementType(type);
			INFO("%d %s elements found.", n_removed_elements, typeName.c_str());
		}
	}
	if (matIDArg.isSet()) {
		const std::vector<unsigned> vec_matID = matIDArg.getValue();
		for (auto matID : vec_matID) {
			const std::size_t n_removed_elements = ex.searchByMaterialID(matID);
			INFO("%d elements with material ID %d found.", n_removed_elements, matID);
		}
	}

	// remove the elements and create a new mesh object.
	MeshLib::Mesh const*const new_mesh = ex.removeMeshElements(mesh->getName());

	// write into a file
	FileIO::Legacy::MeshIO meshIO;
	meshIO.setMesh(new_mesh);
	meshIO.writeToFile(mesh_out.getValue());

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}



