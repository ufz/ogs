/**
 * \file
 * \author Karsten Rink
 * \date   2011-05-10
 * \brief  Implementation of the ElementTreeModel class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementTreeModel.h"

#include "GeoLib/AABB.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshInformation.h"
#include "MeshLib/Node.h"
#include "MeshLib/Vtk/VtkMappedMeshSource.h"
#include "Base/TreeItem.h"

namespace
{
template <typename PropertyType>
QList<QVariant> propertyBounds(PropertyType const& property)
{
    auto const bounds = MeshLib::MeshInformation::getValueBounds(property);
    if (bounds.has_value())
    {
        return {"[" + QString::number(bounds->first) + ",",
                QString::number(bounds->second) + "]"};
    }
    // Makeup the same structure of output as in the valid case above.
    return {"[empty,", "empty]"};
}
}  // namespace

/**
 * Constructor.
 */
ElementTreeModel::ElementTreeModel(QObject* parent) : TreeModel(parent)
{
    QList<QVariant> rootData;
    delete _rootItem;
    rootData << "Name"
             << "Type"
             << ""
             << "";
    _rootItem = new TreeItem(rootData, nullptr);
}

ElementTreeModel::~ElementTreeModel() = default;

void ElementTreeModel::setElement(
    vtkUnstructuredGridAlgorithm const* const grid, const unsigned elem_index)
{
    beginResetModel();

    _mesh_source = grid;
    this->clearView();

    auto const* const source =
        dynamic_cast<MeshLib::VtkMappedMeshSource const* const>(grid);

    if (!source)
    {
        return;
    }

    const MeshLib::Mesh* mesh = source->GetMesh();
    const MeshLib::Element* elem = mesh->getElement(elem_index);

    QList<QVariant> elemData;
    elemData << "Element " + QString::number(elem_index) << ""
             << ""
             << "";
    auto* elemItem = new TreeItem(elemData, _rootItem);
    _rootItem->appendChild(elemItem);

    QList<QVariant> typeData;
    typeData << "Element Type: "
             << QString::fromStdString(
                    MeshElemType2String(elem->getGeomType()));
    auto* typeItem = new TreeItem(typeData, elemItem);
    elemItem->appendChild(typeItem);

    auto const mat_ids = materialIDs(*mesh);
    QString matIdString = !mat_ids ? QString("not defined")
                                   : QString::number((*mat_ids)[elem->getID()]);
    QList<QVariant> materialData;
    materialData << "MaterialID: " << matIdString;
    auto* matItem = new TreeItem(materialData, elemItem);
    elemItem->appendChild(matItem);

    QList<QVariant> volData;
    volData << "Area/Volume: "
            << QString::number(mesh->getElement(elem_index)->getContent());
    auto* volItem = new TreeItem(volData, elemItem);
    elemItem->appendChild(volItem);

    QList<QVariant> nodeListData;
    nodeListData << "Nodes"
                 << ""
                 << ""
                 << "";
    auto* nodeListItem = new TreeItem(nodeListData, elemItem);
    elemItem->appendChild(nodeListItem);

    // const std::vector<MeshLib::Node*> nodes_vec = grid->getNodes();
    std::size_t nElemNodes = elem->getNumberOfBaseNodes();
    for (std::size_t i = 0; i < nElemNodes; i++)
    {
        const MeshLib::Node* node = elem->getNode(i);
        QList<QVariant> nodeData;
        nodeData << "Node " + QString::number(node->getID())
                 << QString::number((*node)[0], 'f', 6)
                 << QString::number((*node)[1], 'f', 6)
                 << QString::number((*node)[2], 'f', 6);
        auto* nodeItem = new TreeItem(nodeData, nodeListItem);
        nodeListItem->appendChild(nodeItem);
    }
    endResetModel();
}

void ElementTreeModel::clearView()
{
    beginResetModel();
    _rootItem->removeChildren(0, _rootItem->childCount());
    endResetModel();
}

void ElementTreeModel::setMesh(MeshLib::Mesh const& mesh)
{
    beginResetModel();

    this->clearView();

    QList<QVariant> mesh_name;
    mesh_name << "Name:" << QString::fromStdString(mesh.getName()) << ""
              << ""
              << "";
    auto* name_item = new TreeItem(mesh_name, _rootItem);
    _rootItem->appendChild(name_item);

    QList<QVariant> nodes_number;
    nodes_number << "#Nodes: " << QString::number(mesh.getNumberOfNodes()) << ""
                 << "";
    auto* nodes_item = new TreeItem(nodes_number, _rootItem);
    _rootItem->appendChild(nodes_item);

    QList<QVariant> elements_number;
    elements_number << "#Elements: "
                    << QString::number(mesh.getNumberOfElements()) << ""
                    << "";
    auto* elements_item = new TreeItem(elements_number, _rootItem);
    _rootItem->appendChild(elements_item);

    auto const& n_element_types =
        MeshLib::MeshInformation::getNumberOfElementTypes(mesh);
    for (auto entry : n_element_types)
    {
        QList<QVariant> number_of_element_types;
        number_of_element_types
            << QString::fromStdString(
                   MeshLib::MeshElemType2String(
                       static_cast<MeshLib::MeshElemType>(entry.first)) +
                   "s:")
            << QString::number(entry.second) << ""
            << "";
        auto* type_item = new TreeItem(number_of_element_types, elements_item);
        elements_item->appendChild(type_item);
    }

    QList<QVariant> bounding_box;
    bounding_box << "Bounding Box"
                 << ""
                 << ""
                 << "";
    auto* aabb_item = new TreeItem(bounding_box, _rootItem);
    _rootItem->appendChild(aabb_item);

    const GeoLib::AABB aabb(MeshLib::MeshInformation::getBoundingBox(mesh));
    auto const& min = aabb.getMinPoint();
    auto const& max = aabb.getMaxPoint();

    QList<QVariant> min_aabb;
    min_aabb << "Min:" << QString::number(min[0], 'f')
             << QString::number(min[1], 'f') << QString::number(min[2], 'f');
    auto* min_item = new TreeItem(min_aabb, aabb_item);
    aabb_item->appendChild(min_item);

    QList<QVariant> max_aabb;
    max_aabb << "Max:" << QString::number(max[0], 'f')
             << QString::number(max[1], 'f') << QString::number(max[2], 'f');
    auto* max_item = new TreeItem(max_aabb, aabb_item);
    aabb_item->appendChild(max_item);

    QList<QVariant> edges;
    edges << "Edge Length: "
          << "[" + QString::number(mesh.getMinEdgeLength(), 'f') + ","
          << QString::number(mesh.getMaxEdgeLength(), 'f') + "]"
          << "";
    auto* edge_item = new TreeItem(edges, _rootItem);
    _rootItem->appendChild(edge_item);

    for (auto [name, property] : mesh.getProperties())
    {
        QList<QVariant> array_info{QString::fromStdString(name) + ": "};

        if (auto p = dynamic_cast<MeshLib::PropertyVector<double>*>(property))
        {
            array_info.append(propertyBounds(*p));
        }
        else if (auto p =
                     dynamic_cast<MeshLib::PropertyVector<float>*>(property))
        {
            array_info.append(propertyBounds(*p));
        }
        else if (auto p = dynamic_cast<MeshLib::PropertyVector<int>*>(property))
        {
            array_info.append(propertyBounds(*p));
        }
        else if (auto p =
                     dynamic_cast<MeshLib::PropertyVector<unsigned>*>(property))
        {
            array_info.append(propertyBounds(*p));
        }
        else if (auto p =
                     dynamic_cast<MeshLib::PropertyVector<long>*>(property))
        {
            array_info.append(propertyBounds(*p));
        }
        else if (auto p = dynamic_cast<MeshLib::PropertyVector<long long>*>(
                     property))
        {
            array_info.append(propertyBounds(*p));
        }
        else if (auto p = dynamic_cast<MeshLib::PropertyVector<unsigned long>*>(
                     property))
        {
            array_info.append(propertyBounds(*p));
        }
        else if (auto p =
                     dynamic_cast<MeshLib::PropertyVector<unsigned long long>*>(
                         property))
        {
            array_info.append(propertyBounds(*p));
        }
        else if (auto p = dynamic_cast<MeshLib::PropertyVector<std::size_t>*>(
                     property))
        {
            array_info.append(propertyBounds(*p));
        }
        else if (auto p =
                     dynamic_cast<MeshLib::PropertyVector<char>*>(property))
        {
            array_info.append(propertyBounds(*p));
        }
        else
        {  // Unhandled property vector type.
            array_info << "[ ?"
                       << "? ]"
                       << "";
        }
        _rootItem->appendChild(new TreeItem(array_info, _rootItem));
    }

    endResetModel();
}
