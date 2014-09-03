/**
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "Elements/Element.h"
#include "MeshEditing/ElementValueModification.h"

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Edit material IDs of mesh elements.", ' ', "0.1");
	TCLAP::SwitchArg replaceArg("r", "replace", "replace material IDs", false);
	TCLAP::SwitchArg condenseArg("c", "condense", "condense material IDs", false);
	cmd.xorAdd(replaceArg, condenseArg);
	TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
	                                     "the name of the file containing the input mesh", true,
	                                     "", "file name");
	cmd.add(mesh_in);
	TCLAP::ValueArg<std::string> mesh_out("o", "mesh-output-file",
	                                      "the name of the file the mesh will be written to", true,
	                                      "", "file name");
	cmd.add(mesh_out);
	TCLAP::MultiArg<unsigned> matIDArg("m", "current-material-id",
	                                      "current material id to be replaced", false, "number");
	cmd.add(matIDArg);
	TCLAP::ValueArg<unsigned> newIDArg("n", "new-material-id",
	                                      "new material id", false, 0, "number");
	cmd.add(newIDArg);
	cmd.parse(argc, argv);

	if (!replaceArg.isSet() && !condenseArg.isSet()) {
		INFO("Please select editing mode: -r or -c");
		return 0;
	} else if (replaceArg.isSet() && condenseArg.isSet()) {
		INFO("Please select only one editing mode: -r or -c");
		return 0;
	} else if (replaceArg.isSet()) {
		if (!matIDArg.isSet() || !newIDArg.isSet()) {
			INFO("current and new material IDs must be provided for relplacement");
			return 0;
		}
	}

	MeshLib::Mesh* mesh (FileIO::readMeshFromFile(mesh_in.getValue()));
	INFO("Mesh read: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());

	if (condenseArg.isSet()) {
		INFO("Condensing material ID...");
		MeshLib::ElementValueModification::condense(*mesh);
	} else if (replaceArg.isSet()) {
		INFO("Replacing material ID...");
		const auto vecOldID = matIDArg.getValue();
		const unsigned newID = newIDArg.getValue();
		for (auto oldID : vecOldID) {
			INFO("%d -> %d", oldID, newID);
			MeshLib::ElementValueModification::replace(*mesh, oldID, newID, true);
		}
	}

	// write into a file
	FileIO::Legacy::MeshIO meshIO;
	meshIO.setMesh(mesh);
	meshIO.writeToFile(mesh_out.getValue());

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}

