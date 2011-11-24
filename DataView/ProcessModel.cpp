/**
 * \file ProcessModel.cpp
 * 18/10/2010 KR Initial implementation
 */

// ** INCLUDES **
#include "ProcessItem.h"
#include "CondObjectListItem.h"
#include "CondItem.h"
#include "ProcessItem.h"
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
	ProcessItem* processParent = this->getProcessParent(c->getProcessType());
	if (processParent == NULL)
	{
		ProcessInfo* pcs = new ProcessInfo(c->getProcessType(), c->getProcessPrimaryVariable(), NULL);
		processParent = this->addProcess(pcs);
	}

	CondObjectListItem* condParent = this->getCondParent(processParent, c->getCondType());
	if (condParent == NULL)
		condParent = this->createCondParent(processParent, c->getCondType() /*, c->getAssociatedGeometryName()*/);

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
		std::cout << "Error in ProcessModel::addConditionItem() - Parent object not found..." << std::endl;
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
		std::cout << "Error in ProcessModel::addConditions() - Specified geometrical object \""
		          << condition->getGeoName() << "\" not found in associated geometry..." 
				  << std::endl;
}

void ProcessModel::addConditions(std::vector<FEMCondition*> &conditions)
{
	for (size_t i = 0; i < conditions.size(); i++)
		this->addCondition(conditions[i]);
}

ProcessItem* ProcessModel::addProcess(ProcessInfo *pcs)
{
	if (this->getProcessParent(pcs->getProcessType()) == NULL)
	{
		this->_project.addProcess(pcs);
		QList<QVariant> processData;
		processData << QVariant(QString::fromStdString(FiniteElement::convertProcessTypeToString(pcs->getProcessType()))) << "";
		ProcessItem* process = new ProcessItem(processData, _rootItem, pcs);
		_rootItem->appendChild(process);
		reset();
		return process;
	}
	else
	{
		std::cout << "Error in ProcessModel::addProcess() - " 
			      << FiniteElement::convertProcessTypeToString(pcs->getProcessType()) 
				  << " already exists." << std::endl;
		return NULL;
	}
}

void ProcessModel::removeFEMConditions(const FiniteElement::ProcessType pcs_type, const FEMCondition::CondType cond_type)
{
	ProcessItem* processParent = this->getProcessParent(pcs_type);
	emit conditionsRemoved(this, pcs_type, cond_type);

	if (cond_type != FEMCondition::UNSPECIFIED)
	{
		CondObjectListItem* condParent = getCondParent(processParent, cond_type);
		removeRows(condParent->row(), 1, index(processParent->row(), 0));
	}
	_project.removeConditions(pcs_type, cond_type);
}

void ProcessModel::removeProcess(const FiniteElement::ProcessType type)
{
	const ProcessItem* processParent = this->getProcessParent(type);
	this->_project.removeProcess(type);
	removeRows(processParent->row(), 1, QModelIndex());
}

void ProcessModel::removeAllProcesses()
{
	int nProcesses = _rootItem->childCount();
	for (int i=0; i<nProcesses; i++)
	{
		ProcessItem* item = static_cast<ProcessItem*>(_rootItem->child(i));
		removeProcess(item->getItem()->getProcessType());
	}
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

ProcessItem* ProcessModel::getProcessParent(const FiniteElement::ProcessType type) const
{
	int nLists = _rootItem->childCount();
	for (int i = 0; i < nLists; i++)
		if (static_cast<ProcessItem*>(_rootItem->child(i))->getItem()->getProcessType() == type)
			return static_cast<ProcessItem*>(_rootItem->child(i));

	return NULL;
}

CondObjectListItem* ProcessModel::getCondParent(TreeItem* parent, const FEMCondition::CondType type)
{
	int nLists = parent->childCount();
	for (int i = 0; i < nLists; i++)
		if (dynamic_cast<CondObjectListItem*>(parent->child(i))->getType() == type)
			return dynamic_cast<CondObjectListItem*>(parent->child(i));
	return NULL;
}

CondObjectListItem* ProcessModel::createCondParent(TreeItem* parent, const FEMCondition::CondType cond_type /*, const std::string &geometry_name*/)
{
	QString condType(QString::fromStdString(FEMCondition::condTypeToString(cond_type)));
	QList<QVariant> condData;
	condData << condType << "";
/*
	const std::vector<GEOLIB::Point*>* pnts = _project.getGEOObjects()->getPointVec(geometry_name);
	if (pnts)
	{
*/
		CondObjectListItem* cond = new CondObjectListItem(condData, parent, cond_type, NULL/*pnts*/);
		parent->appendChild(cond);
		/*emit conditionAdded(this, parent->data(0).toString().toStdString(), type);*/
		return cond;
/*
	}

	return NULL;
*/
}

vtkPolyDataAlgorithm* ProcessModel::vtkSource(const FiniteElement::ProcessType pcs_type, const FEMCondition::CondType cond_type)
{
	ProcessItem* processParent = this->getProcessParent(pcs_type);
	if (processParent)
	{
		CondObjectListItem* condParent = this->getCondParent(processParent, cond_type);
		if (condParent)
			return condParent->vtkSource();
	}
	return NULL;
}
