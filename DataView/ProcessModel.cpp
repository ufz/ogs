/**
 * \file ProcessModel.cpp
 * 18/10/2010 KR Initial implementation
 */

// ** INCLUDES **
#include "ProcessItem.h"
#include "CondObjectListItem.h"
#include "ProcessModel.h"
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
	TreeItem* processParent = this->getProcessParent(QString::fromStdString(convertProcessTypeToString(c->getProcessType())), true);
	CondObjectListItem* condParent = this->getCondParent(processParent, c->getCondType());
	if (condParent == NULL)
		condParent = this->createCondParent(processParent, c->getCondType(), c->getAssociatedGeometryName());

	if (condParent)
	{
		QList<QVariant> condData;
		condData << QString::fromStdString(c->getGeoName()) 
			     << QString::fromStdString(c->getGeoTypeAsString());
		CondItem* condItem = new CondItem(condData, condParent, c);
		condParent->appendChild(condItem);
		// add process information
		//QList<QVariant> pcsData;
		//pcsData << QString::fromStdString(convertProcessTypeToString(c->getProcessType()));
		//TreeItem* pcsInfo = new TreeItem(pcsData, condItem);
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

		//condItem->appendChild(pcsInfo);
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

void ProcessModel::addProcess(ProcessInfo *pcs)
{
	QString pcs_name = QString::fromStdString(FiniteElement::convertProcessTypeToString(pcs->getProcessType()));
	TreeItem* exists = this->getProcessParent(pcs_name, true);
	reset();
}

void ProcessModel::removeFEMConditions(const QString &process_name, FEMCondition::CondType type)
{
	TreeItem* processParent = this->getProcessParent(process_name);
	emit conditionsRemoved(this, process_name.toStdString(), type);

	if (type != FEMCondition::UNSPECIFIED)
	{
		TreeItem* condParent = getCondParent(processParent, type);
		removeRows(condParent->row(), 1, index(processParent->row(), 0));
	}
	_project.removeConditions(process_name.toStdString(), type);
}

void ProcessModel::removeProcess(const QString &process_name)
{
	TreeItem* processParent = this->getProcessParent(process_name);
	removeRows(processParent->row(), 1, QModelIndex());
}

void ProcessModel::removeAllProcesses()
{
	int nProcesses = _rootItem->childCount();
	if (nProcesses > 0)
		removeRows(0, nProcesses, QModelIndex());
}

int ProcessModel::getGEOIndex(const std::string &geo_name,
                                GEOLIB::GEOTYPE type,
                                const std::string &obj_name) const
{
	bool exists(false);
	size_t idx(0);
	if (type == GEOLIB::POINT)
		exists = this->_project.getGEOObjects()->getPointVecObj(geo_name)->getElementIDByName(obj_name, idx);
	else if (type == GEOLIB::POLYLINE)
		exists = this->_project.getGEOObjects()->getPolylineVecObj(geo_name)->getElementIDByName(obj_name,idx);
	else if (type == GEOLIB::SURFACE)
		exists = this->_project.getGEOObjects()->getSurfaceVecObj(geo_name)->getElementIDByName(obj_name,idx);

	if (exists)
		return static_cast<int>(idx);
	return -1;
}

TreeItem* ProcessModel::getProcessParent(const QString &processName, bool create_item)
{
	int nLists = _rootItem->childCount();
	for (int i = 0; i < nLists; i++)
		if (_rootItem->child(i)->data(0).toString().compare(processName) == 0)
			return _rootItem->child(i);

	if (create_item)
	{
		QList<QVariant> processData;
		processData << QVariant(processName) << "";
		TreeItem* process = new TreeItem(processData, _rootItem);
		_rootItem->appendChild(process);
		return process;
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

CondObjectListItem* ProcessModel::createCondParent(TreeItem* parent, FEMCondition::CondType type, const std::string &geometry_name)
{
	QString condType(QString::fromStdString(FEMCondition::condTypeToString(type)));
	QList<QVariant> condData;
	condData << condType << "";

	//std::string geo_name = parent->data(0).toString().toStdString();
	const std::vector<GEOLIB::Point*>* pnts = _project.getGEOObjects()->getPointVec(geometry_name);
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
	TreeItem* processParent = this->getProcessParent(QString::fromStdString(name));
	if (processParent)
	{
		CondObjectListItem* condParent = this->getCondParent(processParent, type);
		if (condParent)
			return condParent->vtkSource();
	}
	return NULL;
}
