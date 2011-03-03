/**
 * \file PolylinesModel.cpp
 * 24/9/2009 LB Initial implementation
 * 05/05/2010 KR 2d graphic functionality removed and various layout changes
 *
 * Implementation of PolylinesModel
 */

#include "LinesModel.h"

#include "VtkPolylinesSource.h"
#include "LineEditDialog.h"
#include "AnalyticalGeometry.h"
#include "OGSError.h"

PolylinesModel::PolylinesModel( QString name, const GEOLIB::PolylineVec* polylineVec, QObject* parent /*= 0*/ )
: Model(name, parent), _polylineVec(polylineVec)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Id" << "x" << "y" << "z";
	_rootItem = new TreeItem(rootData, NULL);
	setData(polylineVec, _rootItem);

	this->constructVTKObject();
}

PolylinesModel::~PolylinesModel()
{
	_vtkSource->Delete();
}

int PolylinesModel::columnCount( const QModelIndex& parent /*= QModelIndex()*/ ) const
{
	Q_UNUSED(parent)

	return 4;
}

void PolylinesModel::constructVTKObject()
{
	_vtkSource = VtkPolylinesSource::New();
	VtkPolylinesSource* source = static_cast<VtkPolylinesSource*>(_vtkSource);
	source->SetName(this->_name + " - Lines");
	source->setPolylines(_polylineVec->getVector());
}

void PolylinesModel::setData(const GEOLIB::PolylineVec* polylineVec, TreeItem* parent)
{
	Q_UNUSED(parent)
	const std::vector<GEOLIB::Polyline*> *lines = polylineVec->getVector();

	int nLines = static_cast<int>(lines->size());
	for (int i=0; i<nLines; i++)
	{
		QList<QVariant> line;
		std::string ply_name("");
		line << "Line " + QString::number(i);
		if (polylineVec->getNameOfElementByID(i, ply_name)) line << QString::fromStdString(ply_name);

		GeoTreeItem* lineItem = new GeoTreeItem(line, _rootItem, (*lines)[i]);
		_rootItem->appendChild(lineItem);

		int nPoints = static_cast<int>((*lines)[i]->getNumberOfPoints());
		for (int j=0; j<nPoints; j++)
		{
			QList<QVariant> pnt;
			pnt << static_cast<int>((*lines)[i]->getPointID(j)) << QString::number((*(*(*lines)[i])[j])[0],'f') << QString::number((*(*(*lines)[i])[j])[1],'f') << QString::number((*(*(*lines)[i])[j])[2],'f');
			TreeItem* child = new TreeItem(pnt, lineItem);
			lineItem->appendChild(child);
		}
	}

	reset();
}

void PolylinesModel::updateData()
{
	clearData();
	setData(_polylineVec, _rootItem);
	//TreeModel::updateData();
	this->_vtkSource->Modified();//>Update();
}

void PolylinesModel::callEditPlyDialog()
{
	LineEditDialog lineEdit(*_polylineVec);
	connect(&lineEdit, SIGNAL(connectPolylines(std::vector<size_t>, bool, bool)), this, SLOT(connectPolylineSegments(std::vector<size_t>, bool, bool)));
	lineEdit.exec();
}

void PolylinesModel::connectPolylineSegments(std::vector<size_t> indexlist, bool closePly, bool triangulatePly)
{
	const std::vector<GEOLIB::Polyline*> *polylines = _polylineVec->getVector();
	std::vector<GEOLIB::Polyline*> ply_list;
	for (size_t i=0; i<indexlist.size(); i++)
		ply_list.push_back( (*polylines)[indexlist[i]] );

	// connect polylines
	GEOLIB::Polyline* new_line = GEOLIB::Polyline::contructPolylineFromSegments(ply_list);

	if (new_line)
	{
		// insert result in a new vector of polylines (because the GEOObjects::appendPolylines()-method requires a vector)
		std::vector<GEOLIB::Polyline*> connected_ply;

		if (closePly)
		{
			new_line = GEOLIB::Polyline::closePolyline(*new_line);
/*
			if (triangulatePly)
			{
				std::list<GEOLIB::Triangle> triangles;
				earClippingTriangulationOfPolygon(new_line, triangles);
			}
*/
		}

		connected_ply.push_back(new_line);

		// send result to whoever needs to know
		std::string tmp_str (this->_polylineVec->getName());
		emit requestAppendPolylines(connected_ply, tmp_str);
	}
	else
		OGSError::box("Error connecting polyines.");
}
