/**
 * \file PntsModel.cpp
 * 24/9/2009 LB Initial implementation
 * 05/05/2010 KR 2d graphic functionality removed and various layout changes
 *
 * Implementation of PntsModel
 */

// ** INCLUDES **
#include "PntsModel.h"
#include "TreeItem.h"
#include "GEOObjects.h"

#include "VtkPointsSource.h"

PntsModel::PntsModel( QString name, const GEOLIB::PointVec* pntVec, QObject* parent /*= 0*/ )
: Model(name, parent), _pntVec(pntVec)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Id" << "x" << "y" << "z";
	_rootItem = new TreeItem(rootData, NULL);
	setData(pntVec, _rootItem);

	this->constructVTKObject();
}

PntsModel::~PntsModel()
{
	_vtkSource->Delete();
}

int PntsModel::columnCount( const QModelIndex &parent /*= QModelIndex()*/ ) const
{
	Q_UNUSED(parent)

	return 4;
}

void PntsModel::constructVTKObject()
{
	_vtkSource = VtkPointsSource::New();
	VtkPointsSource* pntsSource = static_cast<VtkPointsSource*>(_vtkSource);
	pntsSource->SetName(this->_name + " - Points");
	pntsSource->setPoints(_pntVec->getVector());
}

QVariant PntsModel::data( const QModelIndex& index, int role ) const
{
	if (!index.isValid())
		return QVariant();

	if ((size_t)index.row() >= _pntVec->size())
		return QVariant();

	GEOLIB::Point* point = (*_pntVec->getVector())[index.row()];
	if (point == NULL)
		return QVariant();

	switch (role)
	{
	case Qt::DisplayRole:
		switch (index.column())
		{
		case 0:
            return index.row();
		case 1:
			//return (*point)[0];
			return QVariant(QString::number((*point)[0],'f'));
		case 2:
			//return (*point)[1];
			return QVariant(QString::number((*point)[1],'f'));
		case 3:
			//return (*point)[2];
			return QVariant(QString::number((*point)[2],'f'));
		default:
			return QVariant();
		}
		break;

	case Qt::ToolTipRole:
		return QString("(%1, %2, %3)").arg((*point)[0]).arg((*point)[1]).arg((*point)[2]);

	default:
		return QVariant();
	}

}

QVariant PntsModel::headerData( int section, Qt::Orientation orientation, int role /*= Qt::DisplayRole*/ ) const
{
	if (role != Qt::DisplayRole)
		return QVariant();

	if (orientation == Qt::Horizontal)
	{
		switch (section)
		{
		case 0: return "Id";
		case 1: return "x";
		case 2: return "y";
		case 3: return "z";
		default: return QVariant();
		}
	}
	else
		return QString("Row %1").arg(section);
}

void PntsModel::setData(const GEOLIB::PointVec* pntVec, TreeItem* parent)
{
	const std::vector<GEOLIB::Point*> *points = pntVec->getVector();

	int nPoints = static_cast<int>(points->size());
	for (int j=0; j<nPoints; j++)
	{
		QList<QVariant> pnt;
		pnt << j << QString::number((*(*points)[j])[0],'f') << QString::number((*(*points)[j])[1],'f') << QString::number((*(*points)[j])[2],'f');
		TreeItem* child = new TreeItem(pnt, parent);
		parent->appendChild(child);
	}

	reset();
}


bool PntsModel::setData( const QModelIndex& index, const QVariant& value, int role /*= Qt::EditRole*/ )
{

	if (index.isValid() && role == Qt::EditRole)
	{
		GEOLIB::Point* point = (*_pntVec->getVector())[index.row()];
		bool wasConversionSuccesfull = false;
		double x, y, z;
		switch (index.column())
		{
		case 0:
//			id = value.toInt(&wasConversionSuccesfull);
//			if (wasConversionSuccesfull)
//				point->id = id;
//			else
//				return false;
//			break;
		case 1:
			x = value.toDouble(&wasConversionSuccesfull);
			if (wasConversionSuccesfull)
				(*point)[0] = x;
			else
				return false;
			break;
		case 2:
			y = value.toDouble(&wasConversionSuccesfull);
			if (wasConversionSuccesfull)
				(*point)[1] = y;
			else
				return false;
			break;
		case 3:
			z = value.toDouble(&wasConversionSuccesfull);
			if (wasConversionSuccesfull)
				(*point)[2] = z;
			else
				return false;
			break;
		default:
			return false;
		}

		emit dataChanged(index, index);
		return true;
	}

	return false;
}

void PntsModel::updateData()
{
	clearData();
	this->_vtkSource->Delete();

	Model::updateData();
	this->constructVTKObject();
}
