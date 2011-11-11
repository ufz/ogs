/**
 * \file ProcessModel.cpp
 * 18/10/2010 KR Initial implementation
 */

// ** INCLUDES **
#include "CondItem.h"
#include "CondObjectListItem.h"
#include "ConditionModel.h"
#include "FEMCondition.h"
#include "GEOObjects.h"
#include "GeoObject.h"
#include "GeoType.h"

#include <QFileInfo>
#include <vtkPolyDataAlgorithm.h>

ProcessModel::ProcessModel( ProjectData &project, QObject* parent /*= 0*/ )
	: TreeModel(parent), _project(project)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Name" << "Value" << "" << "" << "";
	_rootItem = new TreeItem(rootData, NULL);
}

ProcessModel::~ProcessModel()
{
}

int ProcessModel::columnCount( const QModelIndex &parent /*= QModelIndex()*/ ) const
{
	Q_UNUSED(parent)

	return 2;
}

void ProcessModel::addConditionItem(FEMCondition* c)
{
	TreeItem* geoParent =
	        this->getGEOParent(QString::fromStdString(c->getAssociatedGeometryName()), true);
		CondObjectListItem* condParent = this->getCondParent(geoParent, c->getCondType());
	if (condParent == NULL)
		condParent = this->createCondParent(geoParent, c->getCondType());

	if (condParent)
	{
		QList<QVariant> condData;
		condData << QString::fromStdString(c->getGeoName()) 
			     << QString::fromStdString(c->getGeoTypeAsString());
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
		std::vector<double> dis_value = c->getDisValue();
		TreeItem* disInfo;
		if (!(c->getProcessDistributionType() == FiniteElement::LINEAR ||
		      c->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN))
		{
			for (size_t i = 0; i < dis_value.size(); i++)
				disData << dis_value[i];
			disInfo = new TreeItem(disData, condItem);
		}
		else
		{
			size_t nVals = dis_value.size() / 2;
			disData << static_cast<int>(nVals);
			disInfo = new TreeItem(disData, condItem);
			for (size_t i = 0; i < nVals; i++)
			{
				QList<QVariant> linData;
				linData << dis_value[2 * i] << dis_value[2 * i + 1];
				TreeItem* linInfo = new TreeItem(linData, disInfo);
				disInfo->appendChild(linInfo);
			}
		}

		condItem->appendChild(pcsInfo);
		condItem->appendChild(pvInfo);
		condItem->appendChild(disInfo);

		condParent->addCondition(c);
		reset();
	}
	else
		std::cout << "Error in ConditionModel::addConditionItem() - Parent object not found..." << std::endl;
}

void ProcessModel::addCondition(FEMCondition* condition)
{
	const bool is_domain = (condition->getGeoType() == GEOLIB::GEODOMAIN) ? true : false;

	const GEOLIB::GeoObject* object = condition->getGeoObj();
	if (object == NULL)
	{
		object = _project.getGEOObjects()->getGEOObject(
						 condition->getAssociatedGeometryName(),
						 condition->getGeoType(),
						 condition->getGeoName());
		condition->setGeoObj(object);
	}
	if (object || is_domain)
	{
		_project.addCondition(condition);
		this->addConditionItem(condition);
	}
	else
		std::cout << "Error in ConditionModel::addConditions() - Specified geometrical object "
		          << condition->getGeoName() << " not found in associated geometry..." 
				  << std::endl;
}

void ProcessModel::addConditions(std::vector<FEMCondition*> &conditions)
{
	for (size_t i = 0; i < conditions.size(); i++)
		this->addCondition(conditions[i]);
}
/*
   bool ConditionModel::removeConditionItem(const QModelIndex &idx)
   {
    if (idx.isValid())
    {
        CondItem* item = dynamic_cast<CondItem*>(this->getItem(idx));
        if (item)
        {
            emit conditionRemoved(this, idx);
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
 */

void ProcessModel::removeFEMConditions(const QString &geometry_name, FEMCondition::CondType type)
{
	TreeItem* geoParent = this->getGEOParent(geometry_name);
	emit conditionsRemoved(this, geometry_name.toStdString(), type);

	if ((type == FEMCondition::UNSPECIFIED) || (geoParent->childCount() <= 1)) //remove all conditions for the given geometry
		removeRows(geoParent->row(), 1, QModelIndex());
	else
	{
		TreeItem* condParent = getCondParent(geoParent, type);
		removeRows(condParent->row(), 1, index(geoParent->row(), 0));
	}
	_project.removeConditions(geometry_name.toStdString(), type);
}

int ProcessModel::getGEOIndex(const std::string &geo_name,
                                GEOLIB::GEOTYPE type,
                                const std::string &obj_name) const
{
	bool exists(false);
	size_t idx(0);
	if (type == GEOLIB::POINT)
		exists =
		        this->_project.getGEOObjects()->getPointVecObj(geo_name)->
		        getElementIDByName(
		                obj_name,
		                idx);
	else if (type == GEOLIB::POLYLINE)
		exists =
		        this->_project.getGEOObjects()->getPolylineVecObj(geo_name)->
		        getElementIDByName(
		                obj_name,
		                idx);
	else if (type == GEOLIB::SURFACE)
		exists =
		        this->_project.getGEOObjects()->getSurfaceVecObj(geo_name)->
		        getElementIDByName(
		                obj_name,
		                idx);

	if (exists)
		return idx;
	return -1;
}

TreeItem* ProcessModel::getGEOParent(const QString &geoName, bool create_item)
{
	int nLists = _rootItem->childCount();
	for (int i = 0; i < nLists; i++)
		if (_rootItem->child(i)->data(0).toString().compare(geoName) == 0)
			return _rootItem->child(i);

	if (create_item)
	{
		QList<QVariant> geoData;
		geoData << QVariant(geoName) << "";
		TreeItem* geo = new TreeItem(geoData, _rootItem);
		_rootItem->appendChild(geo);
		return geo;
	}
	return NULL;
}

CondObjectListItem* ProcessModel::getCondParent(TreeItem* parent, FEMCondition::CondType type)
{
	int nLists = parent->childCount();
	for (int i = 0; i < nLists; i++)
		if (dynamic_cast<CondObjectListItem*>(parent->child(i))->getType() == type)
			return dynamic_cast<CondObjectListItem*>(parent->child(i));
	return NULL;
}

CondObjectListItem* ProcessModel::createCondParent(TreeItem* parent, FEMCondition::CondType type)
{
	QString condType(QString::fromStdString(FEMCondition::condTypeToString(type)));
	QList<QVariant> condData;
	condData << condType << "";

	std::string geo_name = parent->data(0).toString().toStdString();
	const std::vector<GEOLIB::Point*>* pnts = _project.getGEOObjects()->getPointVec(geo_name);
	if (pnts)
	{
		CondObjectListItem* cond = new CondObjectListItem(condData, parent, type, pnts);
		parent->appendChild(cond);
		emit conditionAdded(this, parent->data(0).toString().toStdString(), type);
		return cond;
	}
	return NULL;
}

vtkPolyDataAlgorithm* ProcessModel::vtkSource(const std::string &name,
                                                FEMCondition::CondType type)
{
	TreeItem* geoParent = this->getGEOParent(QString::fromStdString(name));
	if (geoParent)
	{
		CondObjectListItem* condParent = this->getCondParent(geoParent, type);
		if (condParent)
			return condParent->vtkSource();
	}
	return NULL;
}
