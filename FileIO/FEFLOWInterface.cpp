/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "FEFLOWInterface.h"

#include <cctype>

#include "logog/include/logog.hpp"
#include <QtXml>

#include "StringTools.h"
#include "FileTools.h"
#include "Polygon.h"
#include "Mesh.h"
#include "Node.h"
#include "Elements/Line.h"
#include "Elements/Hex.h"
#include "Elements/Prism.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Tri.h"

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
			size_t count = 0;
			double dummy_coords[3] = {};
			std::generate(vec_nodes.begin(), vec_nodes.end(), [&]() { return new MeshLib::Node(dummy_coords, count++); });
			line_stream.clear();
		}
		//....................................................................
		// NODE (node index for elements)
		else if (line_string.compare("NODE") == 0)
		{
			assert(!vec_nodes.empty());

			MeshElemType eleType = MeshElemType::INVALID;
			if (fem_dim.n_nodes_of_element == 2)
				eleType = MeshElemType::LINE;
			else if (fem_dim.n_nodes_of_element == 3)
				eleType = MeshElemType::TRIANGLE;
			else if (fem_dim.n_nodes_of_element == 4 && fem_class.dimension == 2)
				eleType = MeshElemType::TRIANGLE;
			else if (fem_dim.n_nodes_of_element == 4 && fem_class.dimension == 3)
				eleType = MeshElemType::TETRAHEDRON;
			else if (fem_dim.n_nodes_of_element == 6 && fem_class.dimension == 3)
				eleType = MeshElemType::PRISM;
			else if (fem_dim.n_nodes_of_element == 8 && fem_class.dimension == 3)
				eleType = MeshElemType::HEXAHEDRON;

			if (eleType == MeshElemType::INVALID) {
				ERR("FEFLOWInterface::readFEFLOWFile(): Unsupported element type with the number of node = %d and dim = %d", fem_dim.n_nodes_of_element, fem_class.dimension);
				std::for_each(vec_nodes.begin(), vec_nodes.end(), [](MeshLib::Node* nod) { delete nod;});
				vec_nodes.clear();
				return nullptr;
			}

			vec_elements.reserve(fem_dim.n_elements);
			for (size_t i=0; i<fem_dim.n_elements; i++)
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
		// SUPERMESH
		else if (line_string.compare("SUPERMESH") == 0)
		{
			readSuperMesh(in, fem_class, &points, &lines);
		}
		//....................................................................
	}
	in.close();

	setMaterialID(fem_class, fem_dim, lines, vec_elements);

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

	std::string project_name(BaseLib::extractBaseNameWithoutExtension(filename));
	if (_geoObjects && points)
		_geoObjects->addPointVec(points, project_name);
	if (_geoObjects && lines)
		_geoObjects->addPolylineVec(lines, project_name);

	return new MeshLib::Mesh(project_name, vec_nodes, vec_elements);
}


void FEFLOWInterface::readNodeCoordinates(std::ifstream &in, const FEM_CLASS &fem_class, const FEM_DIM &fem_dim, std::vector<MeshLib::Node*> &vec_nodes)
{
	const size_t no_nodes_per_layer = (fem_class.dimension == 2) ? fem_dim.n_nodes : fem_dim.n_nodes / (fem_class.n_layers3d + 1);
	const size_t n_lines = no_nodes_per_layer / 12 + 1;
	const size_t n_layers = (fem_class.dimension == 3) ? fem_class.n_layers3d + 1 : 1;
	std::string line_string;
	std::stringstream line_stream;
	double x;
	char dummy_char; // for comma(,)
	// x, y
	for (unsigned k = 0; k < 2; k++)
	{
		// each line
		for (size_t i = 0; i < n_lines; i++)
		{
			getline(in, line_string);
			line_stream.str(line_string);
			for (unsigned j = 0; j < 12; j++)
			{
				if (i * 12 + j >= no_nodes_per_layer)
					break;
				line_stream >> x >> dummy_char;
				for (size_t l = 0; l < n_layers; l++)
				{
					const size_t n = i * 12 + l * no_nodes_per_layer + j;
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

std::vector<size_t> FEFLOWInterface::getNodeList(const std::string &str_ranges)
{
	std::vector<size_t> vec_node_IDs;

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
			const size_t start = vec_node_IDs.back();
			const size_t end = BaseLib::str2number<size_t>(str);
			for (size_t i=start+1; i<end+1; i++)
				vec_node_IDs.push_back(i);
			is_range = false;
		} else {
			BaseLib::trim(str);
			vec_node_IDs.push_back(BaseLib::str2number<size_t>(str));
		}
	}

	return vec_node_IDs;
}

void FEFLOWInterface::readElevation(std::ifstream &in, const FEM_CLASS &fem_class, const FEM_DIM &fem_dim, std::vector<MeshLib::Node*> &vec_nodes)
{
	const size_t no_nodes_per_layer = fem_dim.n_nodes / (fem_class.n_layers3d + 1);
	double z = .0;
	std::string str_nodeList;
	std::string line_string;
	std::stringstream line_stream;
	size_t l = 0;
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
			auto vec_nodeIDs = getNodeList(str_nodeList);
			for (auto n0 : vec_nodeIDs)
			{
				const size_t n = n0 - 1 + l * no_nodes_per_layer;
				(*vec_nodes[n])[2] = z;
			}
			str_nodeList.clear();
		}

		if (mode == 0) {
			break;
		} else if (mode == 1) {
			// slice number
			l++;
			assert(l+1==BaseLib::str2number<size_t>(columns.front()));
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

MeshLib::Element* FEFLOWInterface::readElement(const FEM_DIM &fem_dim, const MeshElemType elem_type, const std::string& line, const std::vector<MeshLib::Node*> &nodes)
{
	std::stringstream ss(line);

	unsigned idx[8];
	for (size_t i = 0; i < fem_dim.n_nodes_of_element; ++i)
		ss >> idx[i];
	MeshLib::Node** ele_nodes = new MeshLib::Node*[fem_dim.n_nodes_of_element];

	switch (elem_type)
	{
		default:
			for (unsigned k(0); k < fem_dim.n_nodes_of_element; ++k)
				ele_nodes[k] = nodes[idx[k]-1];
			break;
		case MeshElemType::HEXAHEDRON:
		case MeshElemType::PRISM:
			const unsigned n_half_nodes = fem_dim.n_nodes_of_element/2;
			for (unsigned k(0); k < n_half_nodes; ++k) {
				ele_nodes[k] = nodes[idx[k+n_half_nodes]-1];
				ele_nodes[k+n_half_nodes] = nodes[idx[k]-1];
			}
			break;
	}

	switch (elem_type)
	{
		case MeshElemType::LINE:
			return new MeshLib::Line(ele_nodes);
		case MeshElemType::TRIANGLE:
			return new MeshLib::Tri(ele_nodes);
		case MeshElemType::QUAD:
			return new MeshLib::Quad(ele_nodes);
		case MeshElemType::TETRAHEDRON:
			return new MeshLib::Tet(ele_nodes);
		case MeshElemType::HEXAHEDRON:
			return new MeshLib::Hex(ele_nodes);
		case MeshElemType::PRISM:
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
	QString str_pt_list1 = xmlEle.text().simplified();
	std::istringstream ss(str_pt_list1.toStdString());
	while (!ss.eof())
	{
		size_t pt_id = 0;
		double pt_xyz[3] = { };
		ss >> pt_id;
		for (int i = 0; i < dim; i++)
			ss >> pt_xyz[i];
		points[pt_id - 1] = new GeoLib::PointWithID(pt_xyz, pt_id);
	}
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
			const size_t n_points = str.toLong();
			QString str_ptId_list = xmlEle.text().simplified();
			{
				GeoLib::Polyline* line = new GeoLib::Polyline(*points);
				lines->push_back(line);
				std::istringstream ss(str_ptId_list.toStdString());
				for (size_t i = 0; i < n_points; i++)
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

void FEFLOWInterface::setMaterialID(const FEM_CLASS &fem_class, const FEM_DIM &fem_dim, const std::vector<GeoLib::Polyline*>* lines, std::vector<MeshLib::Element*> &vec_elements)
{
	if (lines && lines->size() > 1)
	{
		for (size_t i = 0; i < vec_elements.size(); i++)
		{
			MeshLib::Element* e = vec_elements[i];
			const MeshLib::Node gpt = e->getCenterOfGravity();
			size_t matId = 0;
			for (size_t j = 0; j < lines->size(); j++)
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
			e->setValue(matId);
		}
	} else if (fem_class.n_layers3d>0) {
		const size_t no_nodes_per_layer = fem_dim.n_nodes / (fem_class.n_layers3d + 1);
		for (auto* e : vec_elements)
		{
			unsigned e_min_nodeID = std::numeric_limits<unsigned>::max();
			for (size_t i=0; i<e->getNBaseNodes(); i++)
				e_min_nodeID = std::min(e_min_nodeID, e->getNodeIndex(i));
			size_t layer_id = e_min_nodeID / no_nodes_per_layer;
			e->setValue(layer_id);
		}
	}

}

} // end namespace FileIO
