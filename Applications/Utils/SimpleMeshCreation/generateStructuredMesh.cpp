/**
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <algorithm>
#include <string>
#include <vector>

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "BaseLib/LogogSimpleFormatter.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"

// FileIO
#include "FileIO/Legacy/MeshIO.h"

namespace
{
    
/// Get dimension of the mesh element type.
/// @param eleType element type
unsigned getDimension(MeshElemType eleType)
{
    switch (eleType)
    {
        case MeshElemType::LINE:
            return 1;
        case MeshElemType::QUAD:
        case MeshElemType::TRIANGLE:
            return 2;
        case MeshElemType::HEXAHEDRON:
        case MeshElemType::PRISM:
        case MeshElemType::PYRAMID:
        case MeshElemType::TETRAHEDRON:
            return 3;
        default:
            return 0;
    }
}

} // end namespace


int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Structured mesh generator.\n"
			"OpenGeoSys-6 software.\n"
			"Copyright (c) 2012-2014, OpenGeoSys Community "
			"(http://www.opengeosys.org) "
			"Distributed under a Modified BSD License. "
			"See accompanying file LICENSE.txt or "
			"http://www.opengeosys.org/project/license",
		' ',
		BaseLib::BuildInfo::git_version_sha1);

	TCLAP::ValueArg<std::string> mesh_out("o", "mesh-output-file",
	                                      "the name of the file the mesh will be written to", true,
	                                      "", "file name of output mesh");
	cmd.add(mesh_out);
	TCLAP::ValueArg<std::string> eleTypeArg("e", "element-type",
	                                      "element type to be created: line | tri | quad | hex", true, "line", "element type");
	cmd.add(eleTypeArg);
	TCLAP::ValueArg<double> lengthZArg("", "lz",
	                                      "length of a domain in z direction", false, 10.0, "real");
	cmd.add(lengthZArg);
	TCLAP::ValueArg<double> lengthYArg("", "ly",
	                                      "length of a domain in y direction", false, 10.0, "real");
	cmd.add(lengthYArg);
	TCLAP::ValueArg<double> lengthXArg("", "lx",
	                                      "length of a domain in x direction", false, 10.0, "real");
	cmd.add(lengthXArg);
	TCLAP::ValueArg<unsigned> nsubdivZArg("", "nz",
	                                      "the number of subdivision in z direction", false, 10, "integer");
	cmd.add(nsubdivZArg);
	TCLAP::ValueArg<unsigned> nsubdivYArg("", "ny",
	                                      "the number of subdivision in y direction", false, 10, "integer");
	cmd.add(nsubdivYArg);
	TCLAP::ValueArg<unsigned> nsubdivXArg("", "nx",
	                                      "the number of subdivision in x direction", false, 10, "integer");
	cmd.add(nsubdivXArg);

	// parse arguments
	cmd.parse(argc, argv);
	const std::string eleTypeName(eleTypeArg.getValue());
	const MeshElemType eleType = String2MeshElemType(eleTypeName);
	const unsigned dim = getDimension(eleType);

	bool dim_used[3] = {false};
	for (unsigned i=0; i<dim; i++)
		dim_used[i] = true;

	std::vector<TCLAP::ValueArg<double>*> vec_lengthArg = {&lengthXArg, &lengthYArg, &lengthZArg};
	std::vector<TCLAP::ValueArg<unsigned>*> vec_ndivArg = {&nsubdivXArg, &nsubdivYArg, &nsubdivZArg};

	const bool isLengthSet = std::any_of(vec_lengthArg.begin(), vec_lengthArg.end(),
	                                     [&](TCLAP::ValueArg<double> *arg){return arg->isSet();});
	if (!isLengthSet) {
		ERR("Missing input: Length information is not provided at all.");
		return 1;
	} else {
		for (unsigned i=0; i<3; i++) {
			if (dim_used[i] && !vec_lengthArg[i]->isSet()) {
				ERR("Missing input: Length for dimension [%d] is required but missing.", i);
				return 1;
			}
		}
	}

	std::vector<double> length(dim);
	std::vector<unsigned> n_subdivision(dim);
	std::vector<double> vec_dx(dim);
	for (unsigned i=0; i<dim; i++)
	{
		length[i] = vec_lengthArg[i]->getValue();
		n_subdivision[i] = vec_ndivArg[i]->getValue();
		vec_dx[i] = length[i] / n_subdivision[i];
	}

	// generate a mesh
	MeshLib::Mesh* mesh = nullptr;
	switch (eleType)
	{
	case MeshElemType::LINE:
		mesh = MeshLib::MeshGenerator::generateLineMesh(length[0], n_subdivision[0]);
		break;
	case MeshElemType::TRIANGLE:
		mesh = MeshLib::MeshGenerator::generateRegularTriMesh(n_subdivision[0], n_subdivision[1], vec_dx[0], vec_dx[1]);
		break;
	case MeshElemType::QUAD:
		mesh = MeshLib::MeshGenerator::generateRegularQuadMesh(n_subdivision[0], n_subdivision[1], vec_dx[0], vec_dx[1]);
		break;
	case MeshElemType::HEXAHEDRON:
		mesh = MeshLib::MeshGenerator::generateRegularHexMesh(n_subdivision[0], n_subdivision[1], n_subdivision[2], vec_dx[0], vec_dx[1], vec_dx[2]);
		break;
	default:
		ERR("Given element type is not supported.");
		break;
	}

	if (mesh)
	{
		INFO("Mesh created: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());

		// write into a file
		FileIO::Legacy::MeshIO meshIO;
		meshIO.setMesh(mesh);
		meshIO.writeToFile(mesh_out.getValue());

		delete mesh;
	}

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}

