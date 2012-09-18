/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file TetGenInterface.h
 *
 *  Created on 2011-09-12 by Thomas Fischer
 */

#ifndef TETGENINTERFACE_H_
#define TETGENINTERFACE_H_

#include <vector>

// forward declaration of class Node and Element
namespace MeshLib
{
	class Node;
	class Element;
	class Mesh;
}

namespace FileIO
{
/**
 * class TetGenInterface is used to read meshes created by <a href="http://tetgen.berlios.de/">TetGen</a>
 */
class TetGenInterface
{
public:
	TetGenInterface();
	virtual ~TetGenInterface();

	/**
	 * Method reads the TetGen mesh from node file and element file.
	 * @param nodes_fname file name of the nodes file
	 * @param ele_fname file name of the elements file
	 * @return on success the method returns a (pointer to a) CFEMesh, else the method returns NULL
	 */
	MeshLib::Mesh* readTetGenMesh (std::string const& nodes_fname,
	                               std::string const& ele_fname);
	/** in order to have a direct access to the
	 * data structures for nodes and elements we make
	 * class TetGenInterface a friend of the mesh class
	 */
	friend class MeshLib::Mesh;

private:
	/**
	 * Method reads the nodes from stream and stores them in the node vector of the mesh class.
	 * For this purpose it uses methods parseNodesFileHeader() and parseNodes().
	 * @param input the input stream
	 * @return true, if all information is read, false if the method detects an error
	 */
	bool readNodesFromStream(std::ifstream &input);
	/**
	 * Method parses the header of the nodes file created by TetGen
	 * @param line the header is in this string (input)
	 * @param n_nodes number of nodes in the file (output)
	 * @param dim the spatial dimension of the node (output)
	 * @param n_attributes the number of attributes for each node (output)
	 * @param boundary_markers have the nodes boundary information (output)
	 * @return true, if the file header is read, false if the method detects an error
	 */
	bool parseNodesFileHeader(std::string &line, std::size_t& n_nodes, std::size_t& dim,
	                          std::size_t& n_attributes, bool& boundary_markers) const;
	/**
	 * method parses the lines reading the nodes from TetGen nodes file
	 * @param ins the input stream (input)
	 * @param n_nodes the number of nodes to read (input)
	 * @param dim the spatial dimension of the node (input)
	 * @return true, if the nodes are read, false if the method detects an error
	 */
	bool parseNodes(std::ifstream& ins, std::size_t n_nodes, std::size_t dim);

	/**
	 * Method reads the elements from stream and stores them in the element vector of the mesh class.
	 * For this purpose it uses methods parseElementsFileHeader() and parseElements().
	 * @param input the input stream
	 * @return true, if all information is read, false if the method detects an error
	 */
	bool readElementsFromStream(std::ifstream &input);
	/**
	 * Method parses the header of the elements file created by TetGen
	 * @param line
	 * @param n_tets
	 * @param n_nodes_per_tet
	 * @param region_attribute is on output true, if there
	 * @return
	 */
	bool parseElementsFileHeader(std::string &line, std::size_t& n_tets, std::size_t& n_nodes_per_tet,
	                             bool& region_attribute) const;
	/**
	 * Method parses the tetrahedras and put them in the element vector of the mesh class.
	 * @param ins the input stream
	 * @param n_tets the number of tetrahedras that should be read
	 * @param n_nodes_per_tet the number of nodes per tetrahedron
	 * @param region_attribute if region attribute is true, region information is read
	 * @return true, if the tetrahedras are read, false if the method detects an error
	 */
	bool parseElements(std::ifstream& ins, std::size_t n_tets, std::size_t n_nodes_per_tet,
	                   bool region_attribute);

	/**
	 * the nodes later on handed over to the mesh are stored in this vector
	 */
	std::vector<MeshLib::Node*> _nodes;
	/**
	 * the elements (tetrahedrons) later on handed over to the mesh are stored in this vector
	 */
	std::vector<MeshLib::Element*> _elements;
	/**
	 * the value is true if the indexing is zero based, else false
	 */
	bool _zero_based_idx;
};
}

#endif /* TETGENINTERFACE_H_ */
