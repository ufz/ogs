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

std::vector<std::size_t> searchByMaterialID(const std::vector<MeshLib::Element*> & ele_vec, unsigned matID)
{
	std::vector<std::size_t> matchedIDs;
	std::size_t i = 0;
	for (MeshLib::Element* ele : ele_vec) {
		if (ele->getValue()==matID)
			matchedIDs.push_back(i);
		i++;
	}
	return matchedIDs;
}

std::vector<std::size_t> searchByElementType(const std::vector<MeshLib::Element*> & ele_vec, MeshElemType eleType)
{
	std::vector<std::size_t> matchedIDs;
	std::size_t i = 0;
	for (MeshLib::Element* ele : ele_vec) {
		if (ele->getGeomType()==eleType)
			matchedIDs.push_back(i);
		i++;
	}
	return matchedIDs;
}

std::vector<std::size_t> searchByZeroContent(const std::vector<MeshLib::Element*> & ele_vec)
{
	std::vector<std::size_t> matchedIDs;
	std::size_t i = 0;
	for (MeshLib::Element* ele : ele_vec) {
		if (ele->getContent()==.0)
			matchedIDs.push_back(i);
		i++;
	}
	return matchedIDs;
}

void updateUnion(const std::vector<std::size_t> &vec1, std::vector<std::size_t> &vec2)
{
	std::vector<std::size_t> vec_temp(vec1.size() + vec2.size());
	auto it = std::set_union(vec1.begin(), vec1.end(), vec2.begin(), vec2.end(), vec_temp.begin());
	vec_temp.resize(it - vec_temp.begin());
	vec2.assign(vec_temp.begin(), vec_temp.end());
}

std::vector<MeshLib::Element*> excludeElements(const std::vector<MeshLib::Element*> & vec_src_eles, const std::vector<std::size_t> &vec_removed)
{
	std::vector<MeshLib::Element*> vec_dest_eles(vec_src_eles.size() - vec_removed.size());
	std::size_t k=0;
	for (std::size_t i=0; i<vec_src_eles.size(); i++) {
		if (std::find(vec_removed.begin(), vec_removed.end(), i) == vec_removed.end()) {
			vec_dest_eles[k] = vec_src_eles[i];
			k++;
		}
	}
	return vec_dest_eles;
}

void copyNodesElements(	const std::vector<MeshLib::Node*> src_nodes,
						const std::vector<MeshLib::Element*> & src_eles,
						std::vector<MeshLib::Node*> &dst_nodes,
						std::vector<MeshLib::Element*> &dst_eles)
{
	// copy nodes
	dst_nodes.resize(src_nodes.size());
	for (std::size_t i=0; i<dst_nodes.size(); i++) {
		dst_nodes[i] = new MeshLib::Node(*src_nodes[i]);
	}

	// copy elements with new nodes
	dst_eles.resize(src_eles.size());
	for (std::size_t i=0; i<dst_eles.size(); i++) {
		auto* src_ele = src_eles[i];
		auto* dst_ele = src_ele->clone();
		for (unsigned j=0; j<src_ele->getNNodes(); j++) {
			dst_ele->setNode(j, dst_nodes[src_ele->getNode(j)->getID()]);
		}
		dst_eles[i] = dst_ele;
	}
}

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

	MeshLib::Mesh* mesh (FileIO::readMeshFromFile(mesh_in.getValue()));
	INFO("Mesh read: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());

	// search elements IDs to be removed
	std::vector<std::size_t> vec_elementIDs_removed;
	if (zveArg.isSet()) {
		std::vector<std::size_t> vec_matched = searchByZeroContent(mesh->getElements());
		updateUnion(vec_matched, vec_elementIDs_removed);
		INFO("%d zero volume elements found.", vec_matched.size());
	}
	if (eleTypeArg.isSet()) {
		std::vector<std::string> eleTypeNames = eleTypeArg.getValue();
		for (auto typeName : eleTypeNames) {
			MeshElemType type = String2MeshElemType(typeName);
			if (type == MeshElemType::INVALID) continue;
			std::vector<std::size_t> vec_matched = searchByElementType(mesh->getElements(), type);
			updateUnion(vec_matched, vec_elementIDs_removed);
			INFO("%d %s elements found.", vec_matched.size(), typeName.c_str());
		}
	}
	if (matIDArg.isSet()) {
		std::vector<unsigned> vec_matID = matIDArg.getValue();
		for (auto matID : vec_matID) {
			std::vector<std::size_t> vec_matched = searchByMaterialID(mesh->getElements(), matID);
			updateUnion(vec_matched, vec_elementIDs_removed);
			INFO("%d elements with material ID %d found.", vec_matched.size(), matID);
		}
	}

	// remove the elements
	INFO("Removing total %d elements...", vec_elementIDs_removed.size());
	std::vector<MeshLib::Element*> tmp_eles = excludeElements(mesh->getElements(), vec_elementIDs_removed);
	INFO("%d elements remained.", tmp_eles.size());
	std::vector<MeshLib::Node*> new_nodes;
	std::vector<MeshLib::Element*> new_eles;
	copyNodesElements(mesh->getNodes(), tmp_eles, new_nodes, new_eles);

	// create a new mesh object. Unsued nodes are removed while construction
	MeshLib::Mesh* new_mesh(new MeshLib::Mesh(mesh->getName(), new_nodes, new_eles));

	// write into a file
	FileIO::MeshIO meshIO;
	meshIO.setMesh(new_mesh);
	meshIO.writeToFile(mesh_out.getValue());

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}



