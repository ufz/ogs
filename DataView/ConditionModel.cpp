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

void ConditionModel::removeFEMConditions(const std::string &geometry_name, GEOLIB::GEOTYPE type)
{
	for (int j=0; j<this->_rootItem->childCount(); j++)
	{
		TreeItem* parent = _rootItem->child(j);
		for (int i=0; i<parent->childCount(); i++)
		{
			const FEMCondition* cond = static_cast<CondItem*>(parent->child(i))->getItem();
			if (geometry_name.compare(cond->getAssociatedGeometryName()) == 0)
			{
				if (type == GEOLIB::INVALID || type == cond->getGeoType()) 
					parent->removeChildren(i,1);
			}
		}
	}
}