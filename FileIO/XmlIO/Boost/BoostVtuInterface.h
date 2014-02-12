/**
 * \file
 * \author Karsten Rink
 * \date   2012-12-05
 * \brief  Implementation of the BoostVtuInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef BOOSTVTUINTERFACE_H_
#define BOOSTVTUINTERFACE_H_

#include "Writer.h"

#include <string>
#include <vector>

#include "MeshEnums.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/optional.hpp>

class ProjectData;

typedef boost::optional<const boost::property_tree::ptree&> OptionalPtree;

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
class BoostVtuInterface : public Writer
{
public:
	BoostVtuInterface();
	~BoostVtuInterface();

	/// Read an unstructured grid from a VTU file
	static MeshLib::Mesh* readVTUFile(const std::string &file_name);

	/// Decide if the mesh data should be written compressed (default is false).
	void setCompressData(bool flag=true) { _use_compressor = flag; };

	/// Set mesh for writing.
	void setMesh(const MeshLib::Mesh*  mesh) { this->_mesh = const_cast<MeshLib::Mesh*>(mesh); };

protected:
	/// Adds a VTK-DataArray of the given name and datatype to the DOM tree and inserts the data-string at that node
	void addDataArray(boost::property_tree::ptree &parent_node, const std::string &name, const std::string &data_type, const std::string &data, unsigned nComponents = 1);

	bool write(std::ostream& stream);

	std::string _export_name;
	MeshLib::Mesh* _mesh;

private:

	/// Returns the ID used by VTK for a given cell type (e.g. "5" for a triangle, etc.)
	unsigned getVTKElementID(MeshElemType type) const;

	/// Check if the root node really specifies an XML file
	static bool isVTKFile(const boost::property_tree::ptree &vtk_root);

	/// Check if the file really specifies a VTK Unstructured Grid
	static bool isVTKUnstructuredGrid(const boost::property_tree::ptree &vtk_root);

	/// Construct an Element-object from the data given to the method and the data at the current stream position.
	static MeshLib::Element* readElement(std::stringstream &iss, const std::vector<MeshLib::Node*> &nodes, unsigned material, unsigned type);

	static unsigned char* uncompressData(boost::property_tree::ptree const& compressed_data_node);

	/// Get an XML attribute value corresponding to given string from a property tree.
	static const boost::optional<std::string> getXmlAttribute(std::string const& key, boost::property_tree::ptree const& tree);

	/// Find first child of a tree, which is a DataArray and has requested name.
	static const OptionalPtree findDataArray(std::string const& array_name, boost::property_tree::ptree const& tree);

	bool _use_compressor;
};

}

#endif /* BOOSTVTUINTERFACE_H_ */
