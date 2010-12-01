/**
 * \file PolylinesModel.cpp
 * 24/9/2009 LB Initial implementation
 * 05/05/2010 KR 2d graphic functionality removed and various layout changes
 *
 * Implementation of PolylinesModel
 */

#include "LinesModel.h"

#include "VtkPolylinesSource.h"


PolylinesModel::PolylinesModel( QString name, const GEOLIB::PolylineVec* polylineVec, QObject* parent /*= 0*/ )
: Model(name, parent), _polylineVec(polylineVec)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Id" << "x" << "y" << "z";
	_rootItem = new TreeItem(rootData, NULL);
	setData(polylineVec, _rootItem);

	_vtkSource = VtkPolylinesSource::New();
	VtkPolylinesSource* source = static_cast<VtkPolylinesSource*>(_vtkSource);
	source->SetName(name);
	source->setPolylines(polylineVec->getVector());
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

		TreeItem* lineItem = new TreeItem(line, _rootItem);
		_rootItem->appendChild(lineItem);

		int nPoints = static_cast<int>((*lines)[i]->getNumberOfPoints());
		for (int j=0; j<nPoints; j++)
		{
			QList<QVariant> pnt;
			pnt << j << QString::number((*(*(*lines)[i])[j])[0],'f') << QString::number((*(*(*lines)[i])[j])[1],'f') << QString::number((*(*(*lines)[i])[j])[2],'f');
			TreeItem* child = new TreeItem(pnt, lineItem);
			lineItem->appendChild(child);
		}
	}

	reset();
}


void PolylinesModel::updateData()
{
	clearData();
	TreeModel::updateData();
}

