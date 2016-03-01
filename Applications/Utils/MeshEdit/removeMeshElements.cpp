/**
 * @file   removeMeshElements.cpp
 * @author Norihiro Watanabe
 * @date   2013/10/15
 * @brief  Remove mesh elements
 *
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>

// TCLAP
#include "tclap/CmdLine.h"

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "FileIO/readMeshFromFile.h"
#include "FileIO/writeMeshToFile.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"

int main (int argc, char* argv[])
{
	ApplicationsLib::LogogSetup logog_setup;

	TCLAP::CmdLine cmd("Remove mesh elements.", ' ', "0.1");

	// Bounding box params
	TCLAP::ValueArg<double> zLargeArg("", "z-max", "largest allowed extent in z-dimension",
	                                  false, std::numeric_limits<double>::max(), "value");
	cmd.add(zLargeArg);
	TCLAP::ValueArg<double> zSmallArg("", "z-min", "smallest allowed extent in z-dimension",
	                                  false,  -1 * std::numeric_limits<double>::max(), "value");
	cmd.add(zSmallArg);
	TCLAP::ValueArg<double> yLargeArg("", "y-max", "largest allowed extent in y-dimension",
	                                  false, std::numeric_limits<double>::max(), "value");
	cmd.add(yLargeArg);
	TCLAP::ValueArg<double> ySmallArg("", "y-min", "smallest allowed extent in y-dimension",
	                                   false,  -1 * std::numeric_limits<double>::max(), "value");
	cmd.add(ySmallArg);
	TCLAP::ValueArg<double> xLargeArg("", "x-max", "largest allowed extent in x-dimension",
	                                   false, std::numeric_limits<double>::max(), "value");
	cmd.add(xLargeArg);
	TCLAP::ValueArg<double> xSmallArg("", "x-min", "smallest allowed extent in x-dimension",
	                                  false, -1 * std::numeric_limits<double>::max(), "value");
	cmd.add(xSmallArg);

	// Non-bounding-box params
	TCLAP::SwitchArg zveArg("z", "zero-volume", "remove zero volume elements", false);
	cmd.add(zveArg);
	TCLAP::MultiArg<std::string> eleTypeArg("t", "element-type",
	                                      "element type to be removed", false, "element type");
	cmd.add(eleTypeArg);
	TCLAP::MultiArg<unsigned> matIDArg("m", "material-id",
	                                      "material id", false, "material id");
	cmd.add(matIDArg);

	// I/O params
	TCLAP::ValueArg<std::string> mesh_out("o", "mesh-output-file",
	                                      "the name of the file the mesh will be written to", true,
	                                      "", "file name of output mesh");
	cmd.add(mesh_out);
	TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
	                                     "the name of the file containing the input mesh", true,
	                                     "", "file name of input mesh");
	cmd.add(mesh_in);

	cmd.parse(argc, argv);

	std::unique_ptr<MeshLib::Mesh const> mesh(
	    FileIO::readMeshFromFile(mesh_in.getValue()));
	INFO("Mesh read: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());
	MeshLib::ElementSearch searcher(*mesh);

	// search elements IDs to be removed
	if (zveArg.isSet()) {
		const std::size_t n_removed_elements = searcher.searchByContent();
		INFO("%d zero volume elements found.", n_removed_elements);
	}
	if (eleTypeArg.isSet()) {
		const std::vector<std::string> eleTypeNames = eleTypeArg.getValue();
		for (auto typeName : eleTypeNames) {
			const MeshLib::MeshElemType type = MeshLib::String2MeshElemType(typeName);
			if (type == MeshLib::MeshElemType::INVALID) continue;
			const std::size_t n_removed_elements = searcher.searchByElementType(type);
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

	if (xSmallArg.isSet() || xLargeArg.isSet() ||
	    ySmallArg.isSet() || yLargeArg.isSet() ||
	    zSmallArg.isSet() || zLargeArg.isSet())
	{
		bool aabb_error (false);
		if (xSmallArg.getValue() >= xLargeArg.getValue())
		{
		    ERR ("Minimum x-extent larger than maximum x-extent.");
		    aabb_error = true;
		}
		if (ySmallArg.getValue() >= yLargeArg.getValue())
		{
		    ERR ("Minimum y-extent larger than maximum y-extent.");
		    aabb_error = true;
		}
		if (zSmallArg.getValue() >= zLargeArg.getValue())
		{
		    ERR ("Minimum z-extent larger than maximum z-extent.");
		    aabb_error = true;
		}
		if (aabb_error)
		    return EXIT_FAILURE;

		std::array<MathLib::Point3d, 2> extent({{
			MathLib::Point3d(std::array<double,3>{{xSmallArg.getValue(),
				ySmallArg.getValue(), zSmallArg.getValue()}}),
			MathLib::Point3d(std::array<double,3>{{xLargeArg.getValue(),
				yLargeArg.getValue(), zLargeArg.getValue()}})}});
		const std::size_t n_removed_elements = searcher.searchByBoundingBox(
			GeoLib::AABB(extent.begin(), extent.end()));
		INFO("%d elements found.", n_removed_elements);
	}

	// remove the elements and create a new mesh object.
	std::unique_ptr<MeshLib::Mesh const> new_mesh(MeshLib::removeElements(
	    *mesh, searcher.getSearchedElementIDs(), mesh->getName()));

	// write into a file
	FileIO::writeMeshToFile(*new_mesh, mesh_out.getValue());

	return EXIT_SUCCESS;
}



