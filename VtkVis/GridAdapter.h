/**
 * \file GridAdapter.h
 * 24/03/2010 KR Initial implementation
 *
 */


#ifndef GRIDADAPTER_H
#define GRIDADAPTER_H

// ** INCLUDES **
#include "msh_mesh.h"

class vtkImageData; // For conversion from Image to QuadMesh

namespace Mesh_Group
{
	class CFEMesh;
	class CNode;
}

/**
 * \brief Adapter class to convert FEM Mesh to a representation more suited for visualisation purposes
 */
class GridAdapter
{
public:
	/// An element structure consisting of a number of nodes and a MshElemType
	typedef struct
	{
		MshElemType::type	type;
		size_t				material;
		std::vector<size_t> nodes;
	} Element;


	/// Constructor using a FEM-Mesh Object as source
	GridAdapter(const Mesh_Group::CFEMesh* mesh = NULL);

	/// Constructor using a MSH-file as source
	GridAdapter(const std::string &filename);

	~GridAdapter();

	/// Returns the total number of unique material IDs.
	size_t getNumberOfMaterials() const;

	/// Returns the vector of nodes.
	const std::vector<GEOLIB::Point*> *getNodes() const { return _nodes; }

	/// Returns the vector of elements.
	const std::vector<Element*> *getElements() const { return _elems; }

	/// Return a vector of elements for one material group only.
	const std::vector<Element*> *getElements(size_t matID) const;

	/// Returns the grid as a CFEMesh for use in OGS-FEM
	const Mesh_Group::CFEMesh* getCFEMesh() const;

	/// Returns the name of the mesh.
	const std::string getName() const { return _name; };

	/// Sets the name for the mesh.
	void setName(const std::string &name) { _name = name; };

	/// Converts greyscale image to quad mesh
	static Mesh_Group::CFEMesh* convertImgToMesh(vtkImageData* img, const std::pair<double,double> &origin, const double &scalingFactor);

private:
	/// Converts an FEM Mesh to a list of nodes and elements.
	int convertCFEMesh(const Mesh_Group::CFEMesh* mesh);

	/// Reads a MSH file into a list of nodes and elements.
	int readMeshFromFile(const std::string &filename);

	/// Converts a string to a MshElemType
	const MshElemType::type getElementType(const std::string &t);

	/// Converts a GridAdapter into an CFEMesh.
	const Mesh_Group::CFEMesh* toCFEMesh() const;

	static Mesh_Group::CElem* createElement(size_t node1, size_t node2, size_t node3);

	std::string _name;
	std::vector<GEOLIB::Point*> *_nodes;
	std::vector<Element*> *_elems;
	const Mesh_Group::CFEMesh* _mesh;
};

#endif // GRIDADAPTER_H
