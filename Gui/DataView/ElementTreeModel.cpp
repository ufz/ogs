/**
 * \file
 * \author Karsten Rink
 * \date   2011-05-10
 * \brief  Implementation of the ElementTreeModel class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementTreeModel.h"
#include "OGSError.h"
#include "TreeItem.h"
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "AABB.h"

#include "VtkMeshSource.h"

/**
 * Constructor.
 */
ElementTreeModel::ElementTreeModel( QObject* parent )
	: TreeModel(parent)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Name" << "Type" << "" << "";
	_rootItem = new TreeItem(rootData, nullptr);
}

ElementTreeModel::~ElementTreeModel()
{
}

void ElementTreeModel::setElement(vtkUnstructuredGridAlgorithm const*const grid, const unsigned elem_index)
{
	_mesh_source = grid;
	this->clearView();

	VtkMeshSource const*const source = dynamic_cast<VtkMeshSource const*const>(grid);

	if (source)
	{
		const MeshLib::Mesh* mesh = source->GetMesh();
		const MeshLib::Element* elem = mesh->getElement(elem_index);

		QList<QVariant> elemData;
		elemData << "Element " + QString::number(elem_index) << "" << "" << "";
		TreeItem* elemItem = new TreeItem(elemData, _rootItem);
		_rootItem->appendChild(elemItem);

		QList<QVariant> typeData;
		typeData << "Element Type: " << QString::fromStdString(MeshElemType2String(elem->getGeomType()));
		TreeItem* typeItem = new TreeItem(typeData, elemItem);
		elemItem->appendChild(typeItem);

		QList<QVariant> materialData;
		materialData << "MaterialID: " << QString::number(elem->getValue());
		TreeItem* matItem = new TreeItem(materialData, elemItem);
		elemItem->appendChild(matItem);

		QList<QVariant> volData;
		volData << "Area/Volume: " <<
		QString::number(mesh->getElement(elem_index)->getContent());
		TreeItem* volItem = new TreeItem(volData, elemItem);
		elemItem->appendChild(volItem);

		QList<QVariant> nodeListData;
		nodeListData << "Nodes" << "" << "" << "";
		TreeItem* nodeListItem = new TreeItem(nodeListData, elemItem);
		elemItem->appendChild(nodeListItem);

		//const std::vector<MeshLib::Node*> nodes_vec = grid->getNodes();
		size_t nElemNodes = elem->getNNodes();
		for (size_t i = 0; i < nElemNodes; i++)
		{
			const MeshLib::Node* node = elem->getNode(i);
			QList<QVariant> nodeData;
			nodeData << "Node " + QString::number(node->getID()) <<
			QString::number((*node)[0]) << QString::number((*node)[1]) <<
			QString::number((*node)[2]);
			TreeItem* nodeItem = new TreeItem(nodeData, nodeListItem);
			nodeListItem->appendChild(nodeItem);
		}
		reset();
	}
}

void ElementTreeModel::clearView()
{
	_rootItem->removeChildren(0, _rootItem->childCount());
	reset();
}

void ElementTreeModel::setMesh(vtkUnstructuredGridAlgorithm const*const grid)
{
	VtkMeshSource const*const source = dynamic_cast<VtkMeshSource const*const>(grid);
	MeshLib::Mesh const*const mesh = source->GetMesh();
	this->clearView();

	if (mesh)
	{
		QList<QVariant> mesh_name;
		mesh_name << "Name:" << QString::fromStdString(mesh->getName()) << "" << "" << "";
		TreeItem* name_item = new TreeItem(mesh_name, _rootItem);
		_rootItem->appendChild(name_item);

		QList<QVariant> nodes_number;
		nodes_number << "#Nodes: " << QString::number(mesh->getNNodes()) << "" << "";
		TreeItem* nodes_item = new TreeItem(nodes_number, _rootItem);
		_rootItem->appendChild(nodes_item);

		QList<QVariant> elements_number;
		elements_number << "#Elements: " << QString::number(mesh->getNElements()) << "" << "";
		TreeItem* elements_item = new TreeItem(elements_number, _rootItem);
		_rootItem->appendChild(elements_item);

		std::array<unsigned, 7> n_element_types = { 0, 0, 0, 0, 0, 0, 0 };
		std::array<QString, 7> n_element_names = { "Lines:", "Triangles:", "Quads:", "Tetrahedra:", "Hexahedra:", "Pyramids:", "Prisms:" };
		this->getNumberOfElementTypes(mesh, n_element_types);
		for (std::size_t i=0; i<n_element_types.size(); ++i)
		{
			if (n_element_types[i])
			{
				QList<QVariant> elements_number;
				elements_number << n_element_names[i] << QString::number(n_element_types[i]) << "" << "";
				TreeItem* type_item = new TreeItem(elements_number, elements_item);
				elements_item->appendChild(type_item);
			}
		}

		const std::vector<MeshLib::Node*> nodes (mesh->getNodes());
		GeoLib::AABB<MeshLib::Node> aabb(nodes.begin(), nodes.end());
		const MeshLib::Node min = aabb.getMinPoint();
		const MeshLib::Node max = aabb.getMaxPoint();

		QList<QVariant> bounding_box;
		elements_number << "Bounding Box" << "" << "" << "";
		TreeItem* aabb_item = new TreeItem(bounding_box, _rootItem);
		_rootItem->appendChild(aabb_item);

		QList<QVariant> min_aabb;
		min_aabb << "Min:" << QString::number(min[0], 'f') << QString::number(min[1], 'f') << QString::number(min[2], 'f');
		TreeItem* min_item = new TreeItem(min_aabb, aabb_item);
		aabb_item->appendChild(min_item);

		QList<QVariant> max_aabb;
		max_aabb << "Max:" << QString::number(max[0], 'f') << QString::number(max[1], 'f') << QString::number(max[2], 'f');
		TreeItem* max_item = new TreeItem(max_aabb, aabb_item);
		aabb_item->appendChild(max_item);

		std::pair<unsigned, unsigned> mat_bounds (this->getMatBounds(mesh));
		QList<QVariant> materials;
		materials << "MaterialIDs: " << "[" + QString::number(mat_bounds.first) << QString::number(mat_bounds.second) + "]" << "";
		TreeItem* mat_item = new TreeItem(materials, _rootItem);
		_rootItem->appendChild(mat_item);

		reset();
	}
}

std::pair<unsigned, unsigned> ElementTreeModel::getMatBounds(MeshLib::Mesh const*const mesh) const
{
	unsigned min (std::numeric_limits<unsigned>::max());
	unsigned max (0);
	const std::vector<MeshLib::Element*> elements (mesh->getElements());
	for (auto it = elements.begin(); it != elements.end(); ++it)
	{
		const unsigned value = (*it)->getValue();
		if (value > max) max = value;
		if (value < min) min = value;
	}
	return std::pair<unsigned, unsigned>(min, max);
}

void ElementTreeModel::getNumberOfElementTypes(MeshLib::Mesh const*const mesh, std::array<unsigned, 7> &n_element_types) const
{
	const std::vector<MeshLib::Element*> elements (mesh->getElements());
	for (auto it = elements.begin(); it != elements.end(); ++it)
	{
		MeshElemType t = (*it)->getGeomType();
		if (t == MeshElemType::LINE) n_element_types[0]++;
		if (t == MeshElemType::TRIANGLE) n_element_types[1]++;
		if (t == MeshElemType::QUAD) n_element_types[2]++;
		if (t == MeshElemType::TETRAHEDRON) n_element_types[3]++;
		if (t == MeshElemType::HEXAHEDRON) n_element_types[4]++;
		if (t == MeshElemType::PYRAMID) n_element_types[5]++;
		if (t == MeshElemType::PRISM) n_element_types[6]++;
	}
}