/**
 * \file
 * \author Karsten Rink
 * \date   2011-02-07
 * \brief  Implementation of the GeoTreeModel class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <logog/include/logog.hpp>

#include "GeoObjectListItem.h"
#include "GeoTreeItem.h"
#include "GeoTreeModel.h"
#include "OGSError.h"

#include "GeoLib/Triangle.h"

/**
 * Constructor.
 */
GeoTreeModel::GeoTreeModel( QObject* parent )
    : TreeModel(parent)
{
    QList<QVariant> rootData;
    delete _rootItem;
    rootData << "Id" << "x" << "y" << "z" << "name ";
    _rootItem = new GeoTreeItem(rootData, nullptr, nullptr);
}

GeoTreeModel::~GeoTreeModel() = default;

void GeoTreeModel::addPointList(QString geoName, GeoLib::PointVec const& pointVec)
{
    beginResetModel();

    const std::vector<GeoLib::Point*>* points = pointVec.getVector();

    QList<QVariant> geoData;
    geoData << QVariant(geoName) << "" << "" << "" << "";
    auto* geo(new GeoTreeItem(geoData, _rootItem));
    _lists.push_back(geo);
    _rootItem->appendChild(geo);

    QList<QVariant> pointData;
    pointData << "Points" << "" << "" << "" << "";
    auto* pointList =
        new GeoObjectListItem(pointData, geo, points, GeoLib::GEOTYPE::POINT);
    geo->appendChild(pointList);

    std::size_t nPoints = points->size();

    for (std::size_t j = 0; j < nPoints; j++)
    {
        const GeoLib::Point &pnt(*(*points)[j]);
        QList<QVariant> pnt_data;
        pnt_data.reserve(5);
        pnt_data << static_cast<unsigned>(j)
                 << QString::number(pnt[0], 'f')
                 << QString::number(pnt[1], 'f')
                 << QString::number(pnt[2], 'f')
                 << "";
        pointList->appendChild(new GeoTreeItem(pnt_data,
                                               pointList,
                                               static_cast<const GeoLib::Point*>(&pnt)));
    }

    for (auto pnt = pointVec.getNameIDMapBegin(); pnt != pointVec.getNameIDMapEnd(); ++pnt)
        QVariant pnt_data (pointList->child(pnt->second)->setData(4, QString::fromStdString(pnt->first)));

    INFO("Geometry \"%s\" built. %d points added.", geoName.toStdString().c_str(), nPoints);

    endResetModel();
}

void GeoTreeModel::addPolylineList(QString geoName, GeoLib::PolylineVec const& polylineVec)
{
    beginResetModel();

    int nLists = _rootItem->childCount();
    TreeItem* geo(nullptr);
    for (int i = 0; i < nLists; i++)
    {
        if (_rootItem->child(i)->data(0).toString().compare(geoName) == 0)
            geo = _rootItem->child(i);
    }

    if (geo == nullptr)
    {
        ERR("GeoTreeModel::addPolylineList(): No corresponding geometry for \"%s\" found.", geoName.toStdString().c_str());
        return;
    }

    const std::vector<GeoLib::Polyline*>* lines = polylineVec.getVector();

    QList<QVariant> plyData;
    plyData << "Polylines" << "" << "" << "";
    auto* plyList =
        new GeoObjectListItem(plyData, geo, lines, GeoLib::GEOTYPE::POLYLINE);
    geo->appendChild(plyList);
    this->addChildren(plyList, polylineVec, 0, lines->size());

    endResetModel();
}

void GeoTreeModel::appendPolylines(const std::string &name, GeoLib::PolylineVec const& polylineVec)
{
    for (auto& list : _lists)
    {
        if (name == list->data(0).toString().toStdString())
            for (int j = 0; j < list->childCount(); j++)
            {
                auto* parent = static_cast<GeoObjectListItem*>(list->child(j));
                if (GeoLib::GEOTYPE::POLYLINE == parent->getType())
                {
                    beginResetModel();
                    this->addChildren(parent, polylineVec, parent->childCount(),
                                      polylineVec.getVector()->size());
                    endResetModel();
                    parent->vtkSource()->Modified();
                    return;
                }
            }
    }
    OGSError::box("Error adding polyline to geometry.");
}

void GeoTreeModel::addChildren(GeoObjectListItem* plyList,
                               GeoLib::PolylineVec const& polyline_vec,
                               std::size_t start_index,
                               std::size_t end_index)
{
    const std::vector<GeoLib::Polyline*> lines = *(polyline_vec.getVector());

    for (std::size_t i = start_index; i < end_index; i++)
    {
        QList<QVariant> line_data;
        line_data.reserve(4);
        line_data << "Line " + QString::number(i) << "" << "" << "";

        const GeoLib::Polyline &line(*(lines[i]));
        auto* lineItem(new GeoTreeItem(line_data, plyList, &line));
        plyList->appendChild(lineItem);

        auto nPoints = static_cast<int>(lines[i]->getNumberOfPoints());
        for (int j = 0; j < nPoints; j++)
        {
            const GeoLib::Point pnt(*(line.getPoint(j)));
            QList<QVariant> pnt_data;
            pnt_data.reserve(4);
            pnt_data << static_cast<int>(line.getPointID(j))
                     << QString::number(pnt[0], 'f')
                     << QString::number(pnt[1], 'f')
                     << QString::number(pnt[2], 'f');

            lineItem->appendChild(new TreeItem(pnt_data, lineItem));
        }
    }

    for (auto pnt = polyline_vec.getNameIDMapBegin(); pnt != polyline_vec.getNameIDMapEnd(); ++pnt)
        QVariant pnt_data (plyList->child(pnt->second)->setData(1, QString::fromStdString(pnt->first)));

    INFO("%d polylines added.", end_index - start_index);
}

void GeoTreeModel::addSurfaceList(QString geoName, GeoLib::SurfaceVec const& surfaceVec)
{
    beginResetModel();

    int nLists = _rootItem->childCount();
    TreeItem* geo(nullptr);
    for (int i = 0; i < nLists; i++)
    {
        if (_rootItem->child(i)->data(0).toString().compare(geoName) == 0)
            geo = _rootItem->child(i);
    }

    if (geo == nullptr)
    {
        ERR("GeoTreeModel::addSurfaceList(): No corresponding geometry for \"%s\" found.", geoName.toStdString().c_str());
        return;
    }

    const std::vector<GeoLib::Surface*>* surfaces = surfaceVec.getVector();

    QList<QVariant> sfcData;
    sfcData << "Surfaces" << "" << "" << "";
    auto* sfcList =
        new GeoObjectListItem(sfcData, geo, surfaces, GeoLib::GEOTYPE::SURFACE);
    geo->appendChild(sfcList);
    this->addChildren(sfcList, surfaceVec, 0, surfaces->size());

    endResetModel();
}

void GeoTreeModel::appendSurfaces(const std::string &name, GeoLib::SurfaceVec const& surfaceVec)
{
    for (auto& list : _lists)
    {
        if (name == list->data(0).toString().toStdString())
        {
            int nChildren = list->childCount();
            for (int j = 0; j < nChildren; j++)
            {
                auto* parent = static_cast<GeoObjectListItem*>(list->child(j));
                if (GeoLib::GEOTYPE::SURFACE == parent->getType())
                {
                    beginResetModel();
                    this->addChildren(parent, surfaceVec, parent->childCount(),
                                      surfaceVec.getVector()->size());
                    parent->vtkSource()->Modified();
                    endResetModel();
                    return;
                }
            }
        }
    }
    OGSError::box("Error adding surface to geometry.");
}

void GeoTreeModel::addChildren(GeoObjectListItem* sfcList,
                               GeoLib::SurfaceVec const& surface_vec,
                               std::size_t start_index,
                               std::size_t end_index)
{
    const std::vector<GeoLib::Surface*>* surfaces = surface_vec.getVector();

    const std::vector<GeoLib::Point*> &nodesVec(*((*surfaces)[start_index]->getPointVec()));
    for (std::size_t i = start_index; i < end_index; i++)
    {
        QList<QVariant> surface;
        surface.reserve(4);
        surface << "Surface " + QString::number(i) << "" <<
        "" << "";

        const GeoLib::Surface &sfc(*(*surfaces)[i]);
        auto* surfaceItem(new GeoTreeItem(surface, sfcList, &sfc));
        sfcList->appendChild(surfaceItem);

        auto nElems = static_cast<int>((*surfaces)[i]->getNumberOfTriangles());
        for (int j = 0; j < nElems; j++)
        {
            QList<QVariant> elem;
            elem.reserve(4);
            const GeoLib::Triangle &triangle(*sfc[j]);
            elem << j << static_cast<int>(triangle[0])
                 << static_cast<int>(triangle[1])
                 << static_cast<int>(triangle[2]);
            auto* child(new TreeItem(elem, surfaceItem));
            surfaceItem->appendChild(child);

            for (int k = 0; k < 3; k++)
            {
                QList<QVariant> node;
                node.reserve(4);
                const GeoLib::Point &pnt(*(nodesVec[triangle[k]]));
                node << static_cast<int>(triangle[k])
                     << QString::number(pnt[0], 'f')
                     << QString::number(pnt[1], 'f')
                     << QString::number(pnt[2], 'f');
                child->appendChild(new TreeItem(node, child));
            }
        }
    }

    for (auto pnt = surface_vec.getNameIDMapBegin(); pnt != surface_vec.getNameIDMapEnd(); ++pnt)
        QVariant pnt_data (sfcList->child(pnt->second)->setData(1, QString::fromStdString(pnt->first)));

    INFO("%d surfaces added.", end_index - start_index);
}

void GeoTreeModel::renameGeometry(std::string const& old_name,
                                  std::string const& new_name)
{
    for (auto tree_item_entry : _lists)
    {
        if (old_name == tree_item_entry->data(0).toString().toStdString())
        {
            QVariant new_entry(QString::fromStdString(new_name));
            tree_item_entry->setData(0, new_entry);
            break;
        }
    }
    for (auto tree_item_entry : _lists)
    {
        if (new_name == tree_item_entry->data(0).toString().toStdString())
        {
            INFO("Found tree_item_entry with name '%s'.", new_name.c_str());
        }
    }
}

/**
 * Removes the TreeItem with the given name including all its children
 */
void GeoTreeModel::removeGeoList(const std::string& name, GeoLib::GEOTYPE type)
{
    for (std::size_t i = 0; i < _lists.size(); i++)
        if (name == _lists[i]->data(0).toString().toStdString())
        {
            for (int j = 0; j < _lists[i]->childCount(); j++)
                if (type ==
                    static_cast<GeoObjectListItem*>(_lists[i]->child(j))->getType())
                {
                    QModelIndex index = createIndex(j, 0, _lists[i]->child(j));
                    removeRows(0, _lists[i]->child(j)->childCount(), index);
                    removeRows(j, 1, parent(index));
                    break;
                }
            if (_lists[i]->childCount() == 0)
            {
                _lists.erase(_lists.begin() + i);
                removeRows(i, 1, QModelIndex());
            }
        }
}

vtkPolyDataAlgorithm* GeoTreeModel::vtkSource(const std::string &name, GeoLib::GEOTYPE type) const
{
    std::size_t nLists = _lists.size();
    for (std::size_t i = 0; i < nLists; i++)
    {
        if (name == _lists[i]->data(0).toString().toStdString())
            for (int j = 0; j < _lists[i]->childCount(); j++)
            {
                auto* item =
                    dynamic_cast<GeoObjectListItem*>(_lists[i]->child(j));
                if (item->getType() == type)
                    return item->vtkSource();
            }
    }
    return nullptr;
}

void GeoTreeModel::setNameForItem(const std::string &name,
                                  GeoLib::GEOTYPE type,
                                  std::size_t id,
                                  std::string item_name)
{
    std::string geo_type_str("");
    int col_idx(1);

    switch(type)
    {
    case GeoLib::GEOTYPE::POINT:
        geo_type_str = "Points";
        col_idx = 4; // for points the name is at a different position
        break;
    case GeoLib::GEOTYPE::POLYLINE:
        geo_type_str = "Polylines";
        break;
    case GeoLib::GEOTYPE::SURFACE:
        geo_type_str = "Surfaces";
        break;
    default:
        geo_type_str = "";
    }

    auto it = find_if(_lists.begin(), _lists.end(), [&name](GeoTreeItem* geo) {
        return (name == geo->data(0).toString().toStdString());
    });

    for (int i = 0; i < (*it)->childCount(); i++)
    {
        if (geo_type_str == (*it)->child(i)->data(0).toString().toStdString())
        {
            TreeItem* item = (*it)->child(i)->child(id);
            item->setData(col_idx, QString::fromStdString(item_name));
            break;
        }
    }
}
