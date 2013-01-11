/**
 * \file
 * \author Karsten Rink
 * \date   2012-08-30
 * \brief  Definition of the RapidVtuInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef RAPIDVTUINTERFACE_H_
#define RAPIDVTUINTERFACE_H_

#include "Writer.h"

#include <string>
#include <vector>

#include "MshEnums.h"
#include "RapidXML/rapidxml.hpp"

class ProjectData;

namespace MeshLib {
	class Mesh;
	class Node;
	class Element;
}

namespace FileIO
{

/**
 * \brief Reads and writes VtkXMLUnstructuredGrid-files (vtu) to and from OGS data structures.
 *
 * XML parsing is performed using RapidXML. Note, that this library uses string references instead
 * of copies of strings, i.e. strings used in the DOM tree need to be available as long the DOM tree
 * is needed.
 */
class RapidVtuInterface : public Writer
{
public:
	RapidVtuInterface();
	~RapidVtuInterface();

	/// Read an unstructured grid from a VTU file
	static MeshLib::Mesh* readVTUFile(const std::string &file_name);

	/// Decide if the mesh data should be written compressed (default is false).
	void setCompressData(bool flag=true) { _use_compressor = flag; };

	/// Set mesh for writing.
	void setMesh(const MeshLib::Mesh*  mesh) { this->_mesh = const_cast<MeshLib::Mesh*>(mesh); };

protected:
	/// Adds a VTK-DataArray of the given name and datatype to the DOM tree and inserts the data-string at that node
	rapidxml::xml_node<>* addDataArray(const std::string &name, const std::string &data_type, const std::string &data, unsigned nComponents = 1);

	int write(std::ostream& stream);

	std::string _export_name;
	MeshLib::Mesh* _mesh;
	rapidxml::xml_document<> *_doc;

private:
	/// Returns the ID used by VTK for a given cell type (e.g. "5" for a triangle, etc.)
	unsigned getVTKElementID(MshElemType::type type) const;

	/// Check if the root node really specifies an XML file
	static bool isVTKFile(const rapidxml::xml_node<>* node);

	/// Check if the file really specifies a VTK Unstructured Grid
	static bool isVTKUnstructuredGrid(const rapidxml::xml_node<>* node);

	/// Construct an Element-object from the data given to the method and the data at the current stream position.
	static MeshLib::Element* readElement(std::stringstream &iss, const std::vector<MeshLib::Node*> &nodes, unsigned material, unsigned type);

	static unsigned char* uncompressData(const rapidxml::xml_node<>* node);

	bool _use_compressor;
};
}

#endif /* RAPIDVTUINTERFACE_H_ */
