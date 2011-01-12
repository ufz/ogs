/**
 * \file ConditionModel.cpp
 * 18/10/2010 KR Initial implementation
 */

// ** INCLUDES **
#include "ConditionModel.h"
#include "CondItem.h"
#include "GEOObjects.h"
#include "FEMCondition.h"

#include <QFileInfo>

//#include "rf_bc_new.h"
//#include "rf_ic_new.h"
//#include "rf_st_new.h"

//#include "VtkPointsSource.h"

ConditionModel::ConditionModel( ProjectData &project, QObject* parent /*= 0*/ )
: TreeModel(parent), _project(project)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "FEM Conditions";
	_rootItem = new TreeItem(rootData, NULL);

	QList<QVariant> bcData;
	bcData << "Boundary Conditions" << "";
	_bcParent = new TreeItem(bcData, _rootItem);
	_rootItem->appendChild(_bcParent);

	QList<QVariant> icData;
	icData << "Initial Conditions" << "";
	_icParent = new TreeItem(icData, _rootItem);
	_rootItem->appendChild(_icParent);

	QList<QVariant> stData;
	stData << "Source Terms" << "";
	_stParent = new TreeItem(stData, _rootItem);
	_rootItem->appendChild(_stParent);
}

ConditionModel::~ConditionModel()
{
	delete _icParent;
	delete _bcParent;
	delete _stParent;
}

int ConditionModel::columnCount( const QModelIndex &parent /*= QModelIndex()*/ ) const
{
	Q_UNUSED(parent)

	return 2;
}


void ConditionModel::addCondition(FEMCondition* c)
{
	QList<QVariant> condData;
	condData << QString::fromStdString(c->getGeoName()) << QString::fromStdString(c->getGeoTypeAsString());
	TreeItem* parent(NULL);
	size_t row(0);
	if (c->getCondType() == FEMCondition::INITIAL_CONDITION)       { parent = _icParent; row=0; }
	else if (c->getCondType() == FEMCondition::BOUNDARY_CONDITION) { parent = _bcParent; row=1; }
	else if (c->getCondType() == FEMCondition::SOURCE_TERM)		   { parent = _stParent; row=2; }
	CondItem* condItem = new CondItem(condData, parent, c);
	parent->appendChild(condItem);
	// add process information
	QList<QVariant> pcsData;
	pcsData << QString::fromStdString(convertProcessTypeToString(c->getProcessType()));
	TreeItem* pcsInfo = new TreeItem(pcsData, condItem);
	// add information on primary variable
	QList<QVariant> pvData;
	pvData << QString::fromStdString(convertPrimaryVariableToString(c->getProcessPrimaryVariable()));
	TreeItem* pvInfo = new TreeItem(pvData, condItem);
	// add distribution information
	QList<QVariant> disData;
	disData << QString::fromStdString(convertDisTypeToString(c->getProcessDistributionType()));
	TreeItem* disInfo = new TreeItem(disData, condItem);

	condItem->appendChild(pcsInfo);
	condItem->appendChild(pvInfo);
	condItem->appendChild(disInfo);
	//if (stListItem->vtkSource())
	//	stListItem->vtkSource()->SetName(fi.fileName());
	reset();

	emit condAdded(this, this->index(parent->childCount()-1, 0, this->index(row, 0, QModelIndex())));
}


void ConditionModel::addConditions(std::vector<FEMCondition*> &conditions)
{
	for (size_t i=0; i<conditions.size(); i++)
		addCondition(conditions[i]);
}

bool ConditionModel::removeCondition(const QModelIndex &idx)
{
	if (idx.isValid())
	{
		CondItem* item = static_cast<CondItem*>(this->getItem(idx));
		emit condRemoved(this, idx);
		TreeItem* parent = item->parentItem();
		parent->removeChildren(item->row(),1);
		reset();
		return true;
	}

	std::cout << "MshModel::removeMesh() - Specified index does not exist." << std::endl;
	return false;
}

/*
const std::vector<CSourceTerm*> *ConditionModel::getSourceTerms(const QModelIndex &idx) const
{
	if (idx.isValid())
	{
		MshItem* item = static_cast<MshItem*>(this->getItem(idx));
		return item->getGrid();
	}
	std::cout << "MshModel::getMesh() - Specified index does not exist." << std::endl;

	return NULL;
}
*/


/*
QVariant ConditionModel::data( const QModelIndex& index, int role ) const
{
	if (!index.isValid())
		return QVariant();

	if ((size_t)index.row() >= _pntVec->size())
		return QVariant();

	GEOLIB::Point* point = (*_pntVec)[index.row()];
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
*/

//QVariant ConditionModel::headerData( int section, Qt::Orientation orientation, int role /*= Qt::DisplayRole*/ ) const
//{
//	if (role != Qt::DisplayRole)
//		return QVariant();
//
//	if (orientation == Qt::Horizontal)
//	{
//		switch (section)
//		{
//		case 0: return "Id";
//		case 1: return "x";
//		case 2: return "y";
//		case 3: return "z";
//		default: return QVariant();
//		}
//	}
//	else
//		return QString("Row %1").arg(section);
//}

//void ConditionModel::setData(std::vector<CSourceTerm*> *objects, TreeItem* parent)
//{
//	size_t nObjects = objects->size();
//	for (int j=0; j<nObjects; j++)
//	{
//		QList<QVariant> objData;
//		(*objects)[j]->
//		bool pnt = dynamic_cast
//		pnt << j << QString::number((*(*points)[j])[0],'f') << QString::number((*(*points)[j])[1],'f') << QString::number((*(*points)[j])[2],'f');
//		TreeItem* child = new TreeItem(pnt, parent);
//		parent->appendChild(child);
//	}
//
//	reset();
//}


//bool ConditionModel::setData( const QModelIndex& index, const QVariant& value, int role /*= Qt::EditRole*/ )
//{
//
//	if (index.isValid() && role == Qt::EditRole)
//	{
//		GEOLIB::Point* point = (*_pntVec)[index.row()];
//		bool wasConversionSuccesfull = false;
//		double x, y, z;
//		switch (index.column())
//		{
//		case 0:
////			id = value.toInt(&wasConversionSuccesfull);
////			if (wasConversionSuccesfull)
////				point->id = id;
////			else
////				return false;
////			break;
//		case 1:
//			x = value.toDouble(&wasConversionSuccesfull);
//			if (wasConversionSuccesfull)
//				(*point)[0] = x;
//			else
//				return false;
//			break;
//		case 2:
//			y = value.toDouble(&wasConversionSuccesfull);
//			if (wasConversionSuccesfull)
//				(*point)[1] = y;
//			else
//				return false;
//			break;
//		case 3:
//			z = value.toDouble(&wasConversionSuccesfull);
//			if (wasConversionSuccesfull)
//				(*point)[2] = z;
//			else
//				return false;
//			break;
//		default:
//			return false;
//		}
//
//		emit dataChanged(index, index);
//		return true;
//	}
//
//	return false;
//}

/*
void ConditionModel::updateData()
{
	clearData();
	Model::updateData();
}
*/
