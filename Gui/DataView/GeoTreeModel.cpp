/**
 * \file
 * \author Karsten Rink
 * \date   2011-02-07
 * \brief  Implementation of the GeoTreeModel class.
 *
 * \copyright
 * Copyright (c)  2013, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GeoObjectListItem.h"
#include "GeoTreeItem.h"
#include "GeoTreeModel.h"
#include "OGSError.h"

/**
 * Constructor.
 */
GeoTreeModel::GeoTreeModel( QObject* parent )
	: TreeModel(parent)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Id" << "x" << "y" << "z" << "name ";
	_rootItem = new GeoTreeItem(rootData, NULL, NULL);
}

GeoTreeModel::~GeoTreeModel()
{
}

void GeoTreeModel::addPointList(QString geoName, const GeoLib::PointVec* pointVec)
{
	const std::vector<GeoLib::Point*>* points = pointVec->getVector();

	QList<QVariant> geoData;
	geoData << QVariant(geoName) << "" << "" << "" << "";
	GeoTreeItem* geo (new GeoTreeItem(geoData, _rootItem));
	_lists.push_back(geo);
	_rootItem->appendChild(geo);

	QList<QVariant> pointData;
	pointData << "Points" << "" << "" << "" << "";
	GeoObjectListItem* pointList = new GeoObjectListItem(pointData, geo, points, GeoLib::POINT);
	geo->appendChild(pointList);

	size_t nPoints = points->size();

	for (size_t j = 0; j < nPoints; j++)
	{
		const GeoLib::Point &pnt(*(*points)[j]);
		std::string pnt_name("");
		pointVec->getNameOfElementByID(j, pnt_name);
		QList<QVariant> pnt_data;
		pnt_data << static_cast<unsigned>(j)
			     << QString::number(pnt[0], 'f')
			     << QString::number(pnt[1], 'f')
			     << QString::number(pnt[2], 'f')
			     << QString::fromStdString(pnt_name);
		GeoTreeItem* point(new GeoTreeItem(pnt_data, pointList, static_cast<const GeoLib::Point*>(&pnt)));
		pointList->appendChild(point);
	}

	std::cout << "Geometry \"" << geoName.toStdString() << "\" built." << std::endl;
	std::cout << nPoints << " points added." << std::endl;

	reset();
}

void GeoTreeModel::addPolylineList(QString geoName, const GeoLib::PolylineVec* polylineVec)
{
	int nLists = _rootItem->childCount();
	TreeItem* geo(NULL);
	for (int i = 0; i < nLists; i++)
	{
		if (_rootItem->child(i)->data(0).toString().compare(geoName) == 0)
			geo = _rootItem->child(i);
	}

	if (geo == NULL)
	{
		std::cout << "GeoTreeModel::addPolylineList() - Error: No corresponding geometry found..."
				  << std::endl;
		return;
	}

	const std::vector<GeoLib::Polyline*>* lines = polylineVec->getVector();

	QList<QVariant> plyData;
	plyData << "Polylines" << "" << "" << "";
	GeoObjectListItem* plyList = new GeoObjectListItem(plyData, geo, lines, GeoLib::POLYLINE);
	geo->appendChild(plyList);
	this->addChildren(plyList, polylineVec, 0, lines->size());
	reset();
}

void GeoTreeModel::appendPolylines(const std::string &name, const GeoLib::PolylineVec* polylineVec)
{
	for (size_t i = 0; i < _lists.size(); i++)
	{
		if ( name.compare( _lists[i]->data(0).toString().toStdString() ) == 0 )
			for (int j = 0; j < _lists[i]->childCount(); j++)
			{
				GeoObjectListItem* parent =
				        static_cast<GeoObjectListItem*>(_lists[i]->child(j));
				if (GeoLib::POLYLINE == parent->getType())
				{
					this->addChildren(parent, polylineVec,
					                  parent->childCount(),
					                  polylineVec->getVector()->size());
					reset();
					parent->vtkSource()->Modified();
					return;
				}
			}
	}
	OGSError::box("Error adding polyline to geometry.");
}

void GeoTreeModel::addChildren(GeoObjectListItem* plyList,
                               const GeoLib::PolylineVec* polyline_vec,
                               size_t start_index,
                               size_t end_index)
{
	const std::vector<GeoLib::Polyline*> lines = *(polyline_vec->getVector());

	for (size_t i = start_index; i < end_index; i++)
	{
		QList<QVariant> line_data;
		std::string ply_name("");
		if (polyline_vec->getNameOfElementByID(i, ply_name))
			line_data << "Line " + QString::number(i) << QString::fromStdString(ply_name) << "" << "";
		else line_data << "Line " + QString::number(i) << "" << "" << "";

		const GeoLib::Polyline &line(*(lines[i]));
		GeoTreeItem* lineItem(new GeoTreeItem(line_data, plyList, &line));
		plyList->appendChild(lineItem);

		int nPoints = static_cast<int>(lines[i]->getNumberOfPoints());
		for (int j = 0; j < nPoints; j++)
		{
			const GeoLib::Point pnt(*(line[j]));
			QList<QVariant> pnt_data;
			pnt_data << static_cast<int>(line.getPointID(j))
				     << QString::number(pnt[0], 'f')
					 << QString::number(pnt[1], 'f')
					 << QString::number(pnt[2], 'f');
			TreeItem* child(new TreeItem(pnt_data, lineItem));
			lineItem->appendChild(child);
		}
	}
	std::cout << end_index - start_index << " polylines added." << std::endl;
}

void GeoTreeModel::addSurfaceList(QString geoName, const GeoLib::SurfaceVec* surfaceVec)
{
	int nLists = _rootItem->childCount();
	TreeItem* geo(NULL);
	for (int i = 0; i < nLists; i++)
	{
		if (_rootItem->child(i)->data(0).toString().compare(geoName) == 0)
			geo = _rootItem->child(i);
	}

	if (geo == NULL)
	{
		std::cout << "GeoTreeModel::addSurfaceList() - Error: No corresponding geometry found..."
				  << std::endl;
		return;
	}

	const std::vector<GeoLib::Surface*>* surfaces = surfaceVec->getVector();

	QList<QVariant> sfcData;
	sfcData << "Surfaces" << "" << "" << "";
	GeoObjectListItem* sfcList = new GeoObjectListItem(sfcData, geo, surfaces, GeoLib::SURFACE);
	geo->appendChild(sfcList);
	this->addChildren(sfcList, surfaceVec, 0, surfaces->size());

	reset();
}

void GeoTreeModel::appendSurfaces(const std::string &name, GeoLib::SurfaceVec* surfaceVec)
{
	for (size_t i = 0; i < _lists.size(); i++)
	{
		if ( name.compare( _lists[i]->data(0).toString().toStdString() ) == 0 )
		{
			int nChildren = _lists[i]->childCount();
			for (int j = 0; j < nChildren; j++)
			{
				GeoObjectListItem* parent = static_cast<GeoObjectListItem*>(_lists[i]->child(j));
				if (GeoLib::SURFACE == parent->getType())
				{
					this->addChildren(parent, surfaceVec,
					                  parent->childCount(),
					                  surfaceVec->getVector()->size());
					parent->vtkSource()->Modified();
					reset();
					return;
				}
			}
		}
	}
	OGSError::box("Error adding surface to geometry.");
}

void GeoTreeModel::addChildren(GeoObjectListItem* sfcList,
                               const GeoLib::SurfaceVec* surface_vec,
                               size_t start_index,
                               size_t end_index)
{
	const std::vector<GeoLib::Surface*>* surfaces = surface_vec->getVector();

	const std::vector<GeoLib::Point*> &nodesVec(*((*surfaces)[start_index]->getPointVec()));
	for (size_t i = start_index; i < end_index; i++)
	{
		QList<QVariant> surface;
		std::string sfc_name("");
		surface_vec->getNameOfElementByID(i, sfc_name);
		surface << "Surface " + QString::number(i) << QString::fromStdString(sfc_name) << "" << "";

		const GeoLib::Surface &sfc(*(*surfaces)[i]);
		GeoTreeItem* surfaceItem(new GeoTreeItem(surface, sfcList, &sfc));
		sfcList->appendChild(surfaceItem);

		int nElems = static_cast<int>((*surfaces)[i]->getNTriangles());
		for (int j = 0; j < nElems; j++)
		{
			QList<QVariant> elem;
			const GeoLib::Triangle &triangle(*sfc[j]);
			elem << j << static_cast<int>(triangle[0])
				      << static_cast<int>(triangle[1])
					  << static_cast<int>(triangle[2]);
			TreeItem* child(new TreeItem(elem, surfaceItem));
			surfaceItem->appendChild(child);

			for (int k = 0; k < 3; k++)
			{
				QList<QVariant> node;
				const GeoLib::Point &pnt(*(nodesVec[triangle[k]]));
				node << static_cast<int>(triangle[k])
					 << QString::number(pnt[0], 'f')
					 << QString::number(pnt[1], 'f')
					 << QString::number(pnt[2], 'f');
				TreeItem* nchild(new TreeItem(node, child));
				child->appendChild(nchild);
			}
		}
	}
	std::cout << end_index - start_index << " surfaces added." << std::endl;
}

/**
 * Removes the TreeItem with the given name including all its children
 */
void GeoTreeModel::removeGeoList(const std::string &name, GeoLib::GEOTYPE type)
{
	for (size_t i = 0; i < _lists.size(); i++)
		if ( name.compare( _lists[i]->data(0).toString().toStdString() ) == 0 )
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
	size_t nLists = _lists.size();
	for (size_t i = 0; i < nLists; i++)
	{
		if ( name.compare( _lists[i]->data(0).toString().toStdString() ) == 0 )
			for (int j = 0; j < _lists[i]->childCount(); j++)
			{
				GeoObjectListItem* item =
				        dynamic_cast<GeoObjectListItem*>(_lists[i]->child(j));
				if (item->getType() == type)
					return item->vtkSource();
			}
	}
	return NULL;
}

void GeoTreeModel::setNameForItem(const std::string &name, GeoLib::GEOTYPE type, size_t id, std::string item_name)
{
	int type_idx(0);
	int col_idx(1);

	switch(type)
	{
		case GeoLib::POINT:
			type_idx = 0;
			col_idx = 4; // for points the name is at a different position
			break;
		case GeoLib::POLYLINE:
			type_idx = 1;
			break;
		case GeoLib::SURFACE:
			type_idx = 2;
			break;
		case GeoLib::VOLUME:
			type_idx = 3;
			break;
		default:
			type_idx = -1;
	}

	for (size_t i=0; i<_lists.size(); i++)
	{
		if ( name.compare( _lists[i]->data(0).toString().toStdString() ) == 0 )
		{
			TreeItem* object_list = _lists[i]->child(type_idx);
//			for (int j=0; j<object_list->childCount(); j++)
//			{
				TreeItem* item = object_list->child(/*j*/id);
//				if (static_cast<size_t>(item->data(0).toInt()) == id)
//				{
					item->setData(col_idx, QString::fromStdString(item_name));
//					break;
//				}
//			}
			break;
		}
	}
}
