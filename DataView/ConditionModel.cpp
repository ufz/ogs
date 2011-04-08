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
	rootData << "Name" << "Value";;
	_rootItem = new TreeItem(rootData, NULL);
/*
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
*/
}

ConditionModel::~ConditionModel()
{
}

int ConditionModel::columnCount( const QModelIndex &parent /*= QModelIndex()*/ ) const
{
	Q_UNUSED(parent)

	return 2;
}


void ConditionModel::addCondition(FEMCondition* c)
{
	TreeItem* geoParent = getGEOParent(QString::fromStdString(c->getAssociatedGeometryName()));
	std::string g = c->getAssociatedGeometryName();
	TreeItem* condParent = getCondParent(geoParent, c->getCondType());
	std::string cd = condParent->data(0).toString().toStdString();

	QList<QVariant> condData;
	condData << QString::fromStdString(c->getGeoName()) << QString::fromStdString(c->getGeoTypeAsString());
	CondItem* condItem = new CondItem(condData, condParent, c);
	condParent->appendChild(condItem);
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

	emit condAdded(this, this->index(condParent->childCount()-1, 0, this->index(condParent->row(), 0, QModelIndex())));
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
		CondItem* item = dynamic_cast<CondItem*>(this->getItem(idx));
		if (item)
		{
			emit condRemoved(this, idx);
			TreeItem* parent = item->parentItem();
			if (parent->childCount() <=1)
				this->removeFEMConditions(QString::fromStdString(item->getItem()->getAssociatedGeometryName()), item->getItem()->getCondType());
			else
				parent->removeChildren(item->row(),1);
			reset();
			return true;
		}
	}

	std::cout << "ConditionModel::removeCondition() - Specified index does not exist." << std::endl;
	return false;
}

void ConditionModel::removeFEMConditions(const QString &geometry_name, FEMCondition::CondType type = FEMCondition::UNSPECIFIED)
{
	TreeItem* geoParent = getGEOParent(geometry_name);

	if ((type == FEMCondition::UNSPECIFIED) || (geoParent->childCount() <= 1)) //remove all conditions for the given geometry
	{
		_rootItem->removeChildren(geoParent->row(),1);
		return;
	}

	// remove only certain kind of conditions
	TreeItem* condParent = getCondParent(geoParent, type);
	geoParent->removeChildren(condParent->row(), 1);

}



TreeItem* ConditionModel::getGEOParent(const QString &geoName)
{
	int nLists = _rootItem->childCount();
	for (int i=0; i<nLists; i++)
	{
		if (_rootItem->child(i)->data(0).toString().compare(geoName) == 0)
			return _rootItem->child(i);
	}

	QList<QVariant> geoData;
	geoData << QVariant(geoName) << "";
	TreeItem* geo = new TreeItem(geoData, _rootItem);
	//_lists.push_back(geo);
	_rootItem->appendChild(geo);
	return geo;
}

TreeItem* ConditionModel::getCondParent(TreeItem* parent, FEMCondition::CondType type)
{
	QString condType("");
	if (type == FEMCondition::INITIAL_CONDITION) condType = "Initial Conditions";
	else if (type == FEMCondition::BOUNDARY_CONDITION) condType = "Boundary Conditions";
	else if (type == FEMCondition::SOURCE_TERM)	condType = "Source Terms";
	int nLists = parent->childCount();
	for (int i=0; i<nLists; i++)
	{
		if (parent->child(i)->data(0).toString().compare(condType) == 0)
			return parent->child(i);
	}

	QList<QVariant> condData;
	condData << condType << "";
	TreeItem* cond = new TreeItem(condData, parent);
	parent->appendChild(cond);
	return cond;
}
