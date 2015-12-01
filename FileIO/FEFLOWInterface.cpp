/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "FEFLOWInterface.h"

#include <cctype>
#include <QtXml>

#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polygon.h"

#include "MeshLib/Elements/Elements.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace FileIO
{

MeshLib::Mesh* FEFLOWInterface::readFEFLOWFile(const std::string &filename)
{
	std::ifstream in(filename.c_str());
	if (!in)
	{
		ERR("FEFLOWInterface::readFEFLOWFile(): Could not open file %s.", filename.c_str());
		return nullptr;
	}

	FEM_CLASS fem_class;
	FEM_DIM fem_dim;
	std::vector<GeoLib::Point*>* points = NULL;
	std::vector<GeoLib::Polyline*>* lines = NULL;

	bool isXZplane = false;

	std::vector<MeshLib::Node*> vec_nodes;
	std::vector<MeshLib::Element*> vec_elements;

	std::vector<std::vector<std::size_t>> vec_elementsets;

	std::string line_string;
	std::stringstream line_stream;
	while (!in.eof())
	{
		getline(in, line_string);
		//....................................................................
		// CLASS
		if (line_string.find("CLASS") != std::string::npos) // the version number follows afterward, e.g. CLASS (v.5.313)
		{
			getline(in, line_string);
			line_stream.str(line_string);
			// problem class, time mode, problem orientation, dimension, nr. layers for 3D, saturation switch, precision of results, precision of coordinates
			line_stream >> fem_class.problem_class >> fem_class.time_mode >> fem_class.orientation >> fem_class.dimension >> fem_class.n_layers3d;
			line_stream.clear();
		}
		//....................................................................
		// DIMENS
		else if (line_string.compare("DIMENS") == 0)
		{
			// DIMENS
			getline(in, line_string);
			line_stream.str(line_string);
			line_stream >> fem_dim.n_nodes >> fem_dim.n_elements >> fem_dim.n_nodes_of_element >> std::ws;
			// create node pointers with dummy coordinates to create element objects.
			// True coordinates are set later in COOR and ELEV_I.
			vec_nodes.resize(fem_dim.n_nodes);
			std::size_t count = 0;
			double dummy_coords[3] = {};
			std::generate(vec_nodes.begin(), vec_nodes.end(), [&]() { return new MeshLib::Node(dummy_coords, count++); });
			line_stream.clear();
		}
		//....................................................................
		// NODE (node index for elements)
		else if (line_string.compare("NODE") == 0)
		{
			assert(!vec_nodes.empty());

			MeshLib::MeshElemType eleType = MeshLib::MeshElemType::INVALID;
			if (fem_dim.n_nodes_of_element == 2)
				eleType = MeshLib::MeshElemType::LINE;
			else if (fem_dim.n_nodes_of_element == 3)
				eleType = MeshLib::MeshElemType::TRIANGLE;
			else if (fem_dim.n_nodes_of_element == 4 && fem_class.dimension == 2)
				eleType = MeshLib::MeshElemType::QUAD;
			else if (fem_dim.n_nodes_of_element == 4 && fem_class.dimension == 3)
				eleType = MeshLib::MeshElemType::TETRAHEDRON;
			else if (fem_dim.n_nodes_of_element == 6 && fem_class.dimension == 3)
				eleType = MeshLib::MeshElemType::PRISM;
			else if (fem_dim.n_nodes_of_element == 8 && fem_class.dimension == 3)
				eleType = MeshLib::MeshElemType::HEXAHEDRON;

			if (eleType == MeshLib::MeshElemType::INVALID) {
				ERR("FEFLOWInterface::readFEFLOWFile(): Unsupported element type with the number of node = %d and dim = %d", fem_dim.n_nodes_of_element, fem_class.dimension);
				std::for_each(vec_nodes.begin(), vec_nodes.end(), [](MeshLib::Node* nod) { delete nod;});
				vec_nodes.clear();
				return nullptr;
			}

			vec_elements.reserve(fem_dim.n_elements);
			for (std::size_t i=0; i<fem_dim.n_elements; i++)
			{
				getline(in, line_string);
				vec_elements.push_back(readElement(fem_dim, eleType, line_string, vec_nodes));
			}
		}
		//....................................................................
		// COOR
		else if (line_string.compare("COOR") == 0)
		{
			readNodeCoordinates(in, fem_class, fem_dim, vec_nodes);
		}
		//....................................................................
		// ELEV_I
		else if (line_string.compare("ELEV_I") == 0)
		{
			if (fem_class.dimension == 2)
				continue;
			readElevation(in, fem_class, fem_dim, vec_nodes);
		}
		//....................................................................
		// GRAVITY
		else if (line_string.compare("GRAVITY") == 0)
		{
			getline(in, line_string);
			line_stream.str(line_string);
			double vec[3] = { };
			line_stream >> vec[0] >> vec[1] >> vec[2];
			if (vec[0] == 0.0 && vec[1] == -1.0 && vec[2] == 0.0)
				// x-z plane
				isXZplane = true;
			line_stream.clear();
		}
		//....................................................................
		// ELEMENTALSETS
		else if (line_string.compare("ELEMENTALSETS") == 0)
		{
			readELEMENTALSETS(in, vec_elementsets);
		}
		//....................................................................
		// SUPERMESH
		else if (line_string.compare("SUPERMESH") == 0)
		{
			readSuperMesh(in, fem_class, &points, &lines);
		}
		//....................................................................
	}
	in.close();

	std::string project_name(
		BaseLib::extractBaseNameWithoutExtension(filename));
	if (_geoObjects && points)
		_geoObjects->addPointVec(
			std::unique_ptr<std::vector<GeoLib::Point*>>(points), project_name);
	if (_geoObjects && lines)
		_geoObjects->addPolylineVec(
			std::unique_ptr<std::vector<GeoLib::Polyline*>>(lines),
			project_name);

	INFO("setting material IDs");
	auto mesh = new MeshLib::Mesh(project_name, vec_nodes, vec_elements);
	// TODO: Use material id groups
	boost::optional<MeshLib::PropertyVector<int> &> opt_material_ids(
		mesh->getProperties().createNewPropertyVector<int>(
			"MaterialIDs", MeshLib::MeshItemType::Cell, 1)
	);
	if (!opt_material_ids) {
		WARN("Could not create PropertyVector for MaterialIDs in Mesh.");
	}
	else
	{
		opt_material_ids->resize(vec_elements.size());

		if (!vec_elementsets.empty()) {
			for (std::size_t matid =0; matid<vec_elementsets.size(); matid++) {
				auto &eids = vec_elementsets[matid];
				for (auto eid : eids)
					(*opt_material_ids)[eid-1] = matid; // Element IDs given by FEFLOW starts from one!
			}
		} else if (lines && !lines->empty()) {
			for (std::size_t i = 0; i < vec_elements.size(); i++)
			{
				MeshLib::Element* e = vec_elements[i];
				const MeshLib::Node gpt = e->getCenterOfGravity();
				std::size_t matId = 0;
				for (std::size_t j = 0; j < lines->size(); j++)
				{
					GeoLib::Polyline* poly = (*lines)[j];
					if (!poly->isClosed())
						continue;

					GeoLib::Polygon polygon(*poly, true);
					if (polygon.isPntInPolygon(gpt[0], gpt[1], gpt[2]))
					{
						matId = j;
						break;
					}
				}
				(*opt_material_ids)[i] = matId;
			}
		} else if (fem_class.n_layers3d>0) {
			const std::size_t no_nodes_per_layer = fem_dim.n_nodes / (fem_class.n_layers3d + 1);
			for (std::size_t i = 0; i < vec_elements.size(); i++)
			{
				MeshLib::Element* e = vec_elements[i];
				unsigned e_min_nodeID = std::numeric_limits<unsigned>::max();
				for (std::size_t j=0; j<e->getNBaseNodes(); j++)
					e_min_nodeID = std::min(e_min_nodeID, e->getNodeIndex(j));
				std::size_t layer_id = e_min_nodeID / no_nodes_per_layer;
				(*opt_material_ids)[i] = layer_id;
			}
		}
	}

	if (isXZplane)
	{
		for (auto* nod : vec_nodes)
		{
			(*nod)[2] = (*nod)[1];
			(*nod)[1] = 0.0;
		}
		if (points)
		{
			for (auto* pt : *points)
			{
				(*pt)[2] = (*pt)[1];
				(*pt)[1] = .0;
			}
		}
	}

	return mesh;
}


void FEFLOWInterface::readNodeCoordinates(std::ifstream &in, const FEM_CLASS &fem_class, const FEM_DIM &fem_dim, std::vector<MeshLib::Node*> &vec_nodes)
{
	const std::size_t no_nodes_per_layer = (fem_class.dimension == 2) ? fem_dim.n_nodes : fem_dim.n_nodes / (fem_class.n_layers3d + 1);
	assert(no_nodes_per_layer>0);
	const std::size_t n_lines = (no_nodes_per_layer-1) / 12 + 1;
	const std::size_t n_layers = (fem_class.dimension == 3) ? fem_class.n_layers3d + 1 : 1;
	std::string line_string;
	std::stringstream line_stream;
	double x;
	char dummy_char; // for comma(,)
	// x, y
	for (unsigned k = 0; k < 2; k++)
	{
		// each line
		for (std::size_t i = 0; i < n_lines; i++)
		{
			getline(in, line_string);
			line_stream.str(line_string);
			for (unsigned j = 0; j < 12; j++)
			{
				if (i * 12 + j >= no_nodes_per_layer)
					break;
				line_stream >> x >> dummy_char;
				for (std::size_t l = 0; l < n_layers; l++)
				{
					const std::size_t n = i * 12 + l * no_nodes_per_layer + j;
					MeshLib::Node* m_nod = vec_nodes[n];
					if (k == 0)
						(*m_nod)[0] = x;
					else
						(*m_nod)[1] = x;
				}
			}
			line_stream.clear();
		}
	}
}

std::vector<std::size_t> FEFLOWInterface::getIndexList(const std::string &str_ranges)
{
	std::vector<std::size_t> vec_node_IDs;

	// insert space before and after minus for splitting
	std::string str_ranges2(BaseLib::replaceString("-",  " # ", str_ranges));
	BaseLib::trim(str_ranges2);
	auto splitted_str = BaseLib::splitString(str_ranges2, ' ');
	bool is_range = false;
	for (auto str : splitted_str)
	{
		if (str.empty()) continue;
		if (str[0]=='#') {
			is_range = true;
		} else if (is_range) {
			const std::size_t start = vec_node_IDs.back();
			const std::size_t end = BaseLib::str2number<std::size_t>(str);
			for (std::size_t i=start+1; i<end+1; i++)
				vec_node_IDs.push_back(i);
			is_range = false;
		} else {
			BaseLib::trim(str);
			vec_node_IDs.push_back(BaseLib::str2number<std::size_t>(str));
		}
	}

	return vec_node_IDs;
}

void FEFLOWInterface::readElevation(std::ifstream &in, const FEM_CLASS &fem_class, const FEM_DIM &fem_dim, std::vector<MeshLib::Node*> &vec_nodes)
{
	const std::size_t no_nodes_per_layer = fem_dim.n_nodes / (fem_class.n_layers3d + 1);
	double z = .0;
	std::string str_nodeList;
	std::string line_string;
	std::stringstream line_stream;
	std::size_t l = 0;
	unsigned mode = 0; // 0: exit, 1: slice no, 2: elevation value, 3: continued line of mode 2
	int pos_prev_line = 0;
	while (true)
	{
		pos_prev_line = in.tellg();
		getline(in, line_string);

		// check mode
		auto columns = BaseLib::splitString(line_string, ' ');
		if (!in || std::isalpha(line_string[0]))
			mode = 0;
		else if (line_string.empty())
			continue;
		else if (line_string[0]=='\t')
			mode = 3;
		else if (columns.size()==1)
			mode = 1;
		else // columns.size()>1
			mode = 2;

		// process stocked data
		if (mode != 3 && !str_nodeList.empty()) {
			// process previous lines
			auto vec_nodeIDs = getIndexList(str_nodeList);
			for (auto n0 : vec_nodeIDs)
			{
				const std::size_t n = n0 - 1 + l * no_nodes_per_layer;
				(*vec_nodes[n])[2] = z;
			}
			str_nodeList.clear();
		}

		if (mode == 0) {
			break;
		} else if (mode == 1) {
			// slice number
			l++;
			assert(l+1==BaseLib::str2number<std::size_t>(columns.front()));
		} else if (mode == 2) {
			// parse current line
			line_stream.str(line_string);
			line_stream >> z;
			getline(line_stream, str_nodeList);
			BaseLib::trim(str_nodeList, '\t');
			line_stream.clear();
		} else if (mode == 3) {
			// continue reading node range
			BaseLib::trim(line_string, '\t');
			str_nodeList += " " + line_string;
		}
	}

	// move stream position to previous line
	if (std::isalpha(line_string[0]))
		in.seekg(pos_prev_line);
}

MeshLib::Element* FEFLOWInterface::readElement(const FEM_DIM &fem_dim,
	const MeshLib::MeshElemType elem_type, const std::string& line,
	const std::vector<MeshLib::Node*> &nodes)
{
	std::stringstream ss(line);

	unsigned idx[8];
	for (std::size_t i = 0; i < fem_dim.n_nodes_of_element; ++i)
		ss >> idx[i];
	MeshLib::Node** ele_nodes = new MeshLib::Node*[fem_dim.n_nodes_of_element];

	switch (elem_type)
	{
		default:
			for (unsigned k(0); k < fem_dim.n_nodes_of_element; ++k)
				ele_nodes[k] = nodes[idx[k]-1];
			break;
		case MeshLib::MeshElemType::HEXAHEDRON:
		case MeshLib::MeshElemType::PRISM:
			const unsigned n_half_nodes = fem_dim.n_nodes_of_element/2;
			for (unsigned k(0); k < n_half_nodes; ++k) {
				ele_nodes[k] = nodes[idx[k+n_half_nodes]-1];
				ele_nodes[k+n_half_nodes] = nodes[idx[k]-1];
			}
			break;
	}

	switch (elem_type)
	{
		case MeshLib::MeshElemType::LINE:
			return new MeshLib::Line(ele_nodes);
		case MeshLib::MeshElemType::TRIANGLE:
			return new MeshLib::Tri(ele_nodes);
		case MeshLib::MeshElemType::QUAD:
			return new MeshLib::Quad(ele_nodes);
		case MeshLib::MeshElemType::TETRAHEDRON:
			return new MeshLib::Tet(ele_nodes);
		case MeshLib::MeshElemType::HEXAHEDRON:
			return new MeshLib::Hex(ele_nodes);
		case MeshLib::MeshElemType::PRISM:
			return new MeshLib::Prism(ele_nodes);
		default:
			assert(false);
			return nullptr;
	}
}

void FEFLOWInterface::readPoints(QDomElement &nodesEle, const std::string &tag, int dim, std::vector<GeoLib::Point*> &points)
{
	QDomElement xmlEle = nodesEle.firstChildElement(QString::fromStdString(tag));
	if (xmlEle.isNull())
		return;
	QString str_pt_list1 = xmlEle.text();
	std::istringstream ss(str_pt_list1.toStdString());
	std::string line_str;
	while (!ss.eof())
	{
		std::getline(ss, line_str);
		BaseLib::trim(line_str, ' ');
		if (line_str.empty()) continue;
		std::istringstream line_ss(line_str);
		std::size_t pt_id = 0;
		std::array<double,3> pt_xyz;
		line_ss >> pt_id;
		for (int i = 0; i < dim; i++)
			line_ss >> pt_xyz[i];
		points[pt_id - 1] = new GeoLib::Point(pt_xyz, pt_id);
	}
}

void FEFLOWInterface::readELEMENTALSETS(std::ifstream &in, std::vector<std::vector<std::size_t>> &vec_elementsets)
{
	auto compressSpaces = [](std::string const& str) {
		std::stringstream ss(str);
		std::string new_str;
		std::string word;
		while (ss) {
			ss >> word;
			new_str += " " + word;
		}
		return new_str;
	};

	std::string line_string;
	std::string str_idList;
	std::streampos pos_prev_line = 0;
	while (true)
	{
		pos_prev_line = in.tellg();
		getline(in, line_string);

		unsigned mode = 0;
		if (!in)
			mode = 0; // reached the end of the file
		else if (line_string.empty())
			continue; // skip and see what comes next
		else if (std::isalpha(line_string[0]))
			mode = 0; // reached the next section
		else if (line_string[0] == ' ')
			mode = 1; // start of the element set definition
		else if (line_string[0] == '\t')
			mode = 2; // continue the definition
		else
		{
			ERR("Failed during parsing of an ELEMENTALSETS section in a FEFLOW file");
			break;
		}

		if (mode!=2 && !str_idList.empty()) {
			vec_elementsets.push_back(getIndexList(str_idList));
			str_idList.clear();
		}

		if (mode == 0) {
			break;
		} else if (mode == 1) {
			// starting a new set
			std::string set_name;
			std::string ids;
			BaseLib::trim(line_string, ' ');
			if (line_string[0]=='"') { // multiple words
				auto pos = line_string.find_last_of('"');
				set_name = line_string.substr(1, pos-1); // without quotation
				ids = line_string.substr(pos+1);
			} else { // single word
				auto pos = line_string.find_first_of(' ');
				set_name = line_string.substr(0, pos);
				ids = line_string.substr(pos+1);
			}
			INFO("Found an element group - %s", set_name.data());
			str_idList += compressSpaces(ids);
		} else {
			// continue reading a element ids
			BaseLib::trim(line_string, '\t');
			str_idList += compressSpaces(line_string);
		}
	}
	// move stream position to previous line
	if (std::isalpha(line_string[0]))
		in.seekg(pos_prev_line);

}

//
void FEFLOWInterface::readSuperMesh(std::ifstream &in, const FEM_CLASS &fem_class, std::vector<GeoLib::Point*>** p_points, std::vector<GeoLib::Polyline*>** p_lines)
{
	// get XML strings
	std::ostringstream oss;
	std::string line_string;
	while (true)
	{
		getline(in, line_string);
		BaseLib::trim(line_string);
		oss << line_string << "\n";
		if (line_string.find("</supermesh>") != std::string::npos)
			break;
	}
	const QString strXML(oss.str().c_str());

	// convert string to XML
	QDomDocument doc;
	if (!doc.setContent(strXML))
	{
		ERR("FEFLOWInterface::readSuperMesh(): Illegal XML format error");
		return;
	}

	// get geometry data from XML
	QDomElement docElem = doc.documentElement(); // #supermesh
	// #nodes
	*p_points = new std::vector<GeoLib::Point*>();
	std::vector<GeoLib::Point*>* points = *p_points;
	QDomElement nodesEle = docElem.firstChildElement("nodes");
	if (nodesEle.isNull())
		return;

	{
		const QString str = nodesEle.attribute("count");
		const long n_points = str.toLong();
		points->resize(n_points);
		//fixed
		readPoints(nodesEle, "fixed", fem_class.dimension, *points);
		readPoints(nodesEle, "linear", fem_class.dimension, *points);
		readPoints(nodesEle, "parabolic", fem_class.dimension, *points);
	}

	// #polygons
	*p_lines = new std::vector<GeoLib::Polyline*>();
	std::vector<GeoLib::Polyline*>* lines = *p_lines;
	QDomElement polygonsEle = docElem.firstChildElement("polygons");
	if (polygonsEle.isNull())
		return;

	{
		QDomNode child = polygonsEle.firstChild();
		while (!child.isNull())
		{
			if (child.nodeName() != "polygon")
			{
				child = child.nextSibling();
				continue;
			}
			QDomElement xmlEle = child.firstChildElement("nodes");
			if (xmlEle.isNull())
				continue;
			const QString str = xmlEle.attribute("count");
			const std::size_t n_points = str.toLong();
			QString str_ptId_list = xmlEle.text().simplified();
			{
				GeoLib::Polyline* line = new GeoLib::Polyline(*points);
				lines->push_back(line);
				std::istringstream ss(str_ptId_list.toStdString());
				for (std::size_t i = 0; i < n_points; i++)
				{
					int pt_id = 0;
					ss >> pt_id;
					line->addPoint(pt_id - 1);
				}
				line->addPoint(line->getPointID(0));
			}
			child = child.nextSibling();
		}
	}
}

} // end namespace FileIO
