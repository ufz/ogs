/**
 * \file
 * \author Karsten Rink
 * \date   2010-10-18
 * \brief  Implementation of the ProcessModel class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessModel.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// ** INCLUDES **
#include "ProcessItem.h"
#include "CondObjectListItem.h"
#include "CondItem.h"
#include "ProcessItem.h"
#include "FEMCondition.h"
#include "GEOObjects.h"
#include "GeoObject.h"
#include "GeoType.h"
#include "FEMEnums.h"

#include <QFileInfo>
#include <vtkPolyDataAlgorithm.h>

ProcessModel::ProcessModel(ProjectData &project, QObject* parent /*= 0*/) :
		TreeModel(parent), _project(project)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Name" << "Value" << "" << "" << "";
	_rootItem = new TreeItem(rootData, nullptr);
}

ProcessModel::~ProcessModel()
{
}

int ProcessModel::columnCount(const QModelIndex &parent /*= QModelIndex()*/) const
{
	Q_UNUSED(parent)

	return 2;
}

void ProcessModel::addConditionItem(FEMCondition* c)
{
	ProcessItem* processParent = this->getProcessParent(c->getProcessType());
	if (processParent == nullptr)
	{
		ProcessInfo* pcs = new ProcessInfo(c->getProcessType(),
				c->getProcessPrimaryVariable()/* TODO6, nullptr*/);
		processParent = this->addProcess(pcs);
	}

	CondObjectListItem* condParent = this->getCondParent(processParent, c->getCondType());
	if (condParent == nullptr)
		condParent = this->createCondParent(processParent, c);
	else
		condParent->addCondition(c);

	QList<QVariant> condData;
	condData << QString::fromStdString(c->getGeoName())
			<< QString::fromStdString(c->getGeomTypeAsString());
	CondItem* condItem = new CondItem(condData, condParent, c);
	condParent->appendChild(condItem);
	// add information on primary variable
	QList<QVariant> pvData;
	pvData << QString::fromStdString(convertPrimaryVariableToString(c->getProcessPrimaryVariable()));
	TreeItem* pvInfo = new TreeItem(pvData, condItem);
	// add distribution information
	QList<QVariant> disData;
	disData << QString::fromStdString(convertDisTypeToString(c->getProcessDistributionType()));
	std::vector<size_t> dis_nodes = c->getDisNodes();
	std::vector<double> dis_values = c->getDisValues();
	TreeItem* disInfo;
	if (c->getProcessDistributionType() == FiniteElement::CONSTANT ||
	    c->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN ||
		c->getProcessDistributionType() == FiniteElement::NODESCONSTANT)
	{
		disData << dis_values[0];
		disInfo = new TreeItem(disData, condItem);
	}
	else
	{
		size_t nVals = dis_values.size();
		disData << static_cast<int>(nVals);
		disInfo = new TreeItem(disData, condItem);
		for (size_t i = 0; i < nVals; i++)
		{
			QList<QVariant> linData;
			linData << static_cast<int>(dis_nodes[i]) << dis_values[i];
			TreeItem* linInfo = new TreeItem(linData, disInfo);
			disInfo->appendChild(linInfo);
		}
	}

	//condItem->appendChild(pcsInfo);
	condItem->appendChild(pvInfo);
	condItem->appendChild(disInfo);

	reset();
}

void ProcessModel::addCondition(FEMCondition* condition)
{
	// HACK: direct source terms are not domain conditions but they also don't contain geo-object-names
	bool is_domain (false);
	if (condition->getProcessDistributionType() == FiniteElement::NODESCONSTANT ||
		condition->getProcessDistributionType() == FiniteElement::DIRECT)
		is_domain = true;		

	const GeoLib::GeoObject* object = condition->getGeoObj();
	if (object == nullptr)
	{
		object = _project.getGEOObjects()->getGeoObject(condition->getAssociatedGeometryName(),
				condition->getGeomType(), condition->getGeoName());
		condition->setGeoObj(object);
	}
	if (object || is_domain)
	{
		this->addConditionItem(condition);
	}
	else
		ERR("ProcessModel::addConditions(): Specified geometrical object \"%s\" not found in associated geometry.",
				condition->getGeoName().c_str());
}

void ProcessModel::addConditions(std::vector<FEMCondition*> &conditions)
{
	for (size_t i = 0; i < conditions.size(); i++)
		this->addCondition(conditions[i]);
}

ProcessItem* ProcessModel::addProcess(ProcessInfo *pcs)
{
	if (this->getProcessParent(pcs->getProcessType()) == nullptr)
	{
		this->_project.addProcess(pcs);
		QList<QVariant> processData;
		processData
				<< QVariant(
						QString::fromStdString(
								FiniteElement::convertProcessTypeToString(pcs->getProcessType())))
				<< "";
		ProcessItem* process = new ProcessItem(processData, _rootItem, pcs);
		_rootItem->appendChild(process);
		reset();
		return process;
	}
	else
	{
		WARN("ProcessModel::addProcess(): %s  already exists.",
				FiniteElement::convertProcessTypeToString(pcs->getProcessType()).c_str());
		return nullptr;
	}
}

const FEMCondition* ProcessModel::getCondition(const FiniteElement::ProcessType pcs_type, const std::string &geo_name, const FEMCondition::CondType cond_type, const GeoLib::GEOTYPE geo_type, const std::string &cond_name) const
{
	Q_UNUSED(geo_name);
	ProcessItem const*const pcs_item (this->getProcessParent(pcs_type));
	if (!pcs_item)
		return nullptr;
	CondObjectListItem const*const cnd_list = getCondParent(pcs_item, cond_type);
	if (!cnd_list)
		return nullptr;
	for (int i = 0; i < cnd_list->childCount(); i++)
	{
		CondItem* item = static_cast<CondItem*>(cnd_list->child(i));
		if ((item->data(0).toString().toStdString().compare(cond_name) == 0) &&
			(item->data(1).toString().toStdString().compare(GeoLib::convertGeoTypeToString(geo_type)) == 0))
			return item->getItem();
	}
	return nullptr;
}

void ProcessModel::removeConditions(const FiniteElement::ProcessType pcs_type,
		const std::string &geometry_name, const FEMCondition::CondType cond_type)
{
	_project.removeConditions(pcs_type, geometry_name, cond_type);

	while (_rootItem->childCount() > 0)
	{
		ProcessItem* pcs = static_cast<ProcessItem*>(_rootItem->child(0));
		for (int j = 0; j < pcs->childCount(); j++)
			emit conditionsRemoved(this, pcs->getItem()->getProcessType(),
					(static_cast<CondObjectListItem*>(pcs->child(j)))->getType());

		_rootItem->removeChildren(0, 1);
	}

	const std::vector<FEMCondition*> conds = _project.getConditions(FiniteElement::INVALID_PROCESS,
			"", FEMCondition::UNSPECIFIED);
	if (!conds.empty())
	{
		size_t nConds(conds.size());
		for (size_t i = 0; i < nConds; i++)
			this->addConditionItem(conds[i]);
	}
	reset();
}

void ProcessModel::removeProcess(const FiniteElement::ProcessType type)
{
	this->removeConditions(type, "", FEMCondition::UNSPECIFIED);

	const ProcessItem* processParent = this->getProcessParent(type);
	if (processParent)
	{
		this->_project.removeProcess(type);
		removeRows(processParent->row(), 1, QModelIndex());
	}
	reset();
}

void ProcessModel::removeAllProcesses()
{
	int nProcesses = _rootItem->childCount();
	for (int i = 0; i < nProcesses; i++)
	{
		ProcessItem* item = static_cast<ProcessItem*>(_rootItem->child(i));
		removeProcess(item->getItem()->getProcessType());
	}
}

int ProcessModel::getGEOIndex(const std::string &geo_name, GeoLib::GEOTYPE type,
		const std::string &obj_name) const
{
	bool exists(false);
	size_t idx(0);
	if (type == GeoLib::GEOTYPE::POINT)
		exists = this->_project.getGEOObjects()->getPointVecObj(geo_name)->getElementIDByName(
				obj_name, idx);
	else if (type == GeoLib::GEOTYPE::POLYLINE)
		exists = this->_project.getGEOObjects()->getPolylineVecObj(geo_name)->getElementIDByName(
				obj_name, idx);
	else if (type == GeoLib::GEOTYPE::SURFACE)
		exists = this->_project.getGEOObjects()->getSurfaceVecObj(geo_name)->getElementIDByName(
				obj_name, idx);

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
	return nullptr;
}

CondObjectListItem* ProcessModel::getCondParent(TreeItem const*const parent, const FEMCondition::CondType type) const
{
	int nLists = parent->childCount();
	for (int i = 0; i < nLists; i++)
		if (dynamic_cast<CondObjectListItem*>(parent->child(i))->getType() == type)
			return dynamic_cast<CondObjectListItem*>(parent->child(i));
	return nullptr;
}

CondObjectListItem* ProcessModel::createCondParent(ProcessItem* parent, FEMCondition* cond)
{
	QString condType(QString::fromStdString(FEMCondition::condTypeToString(cond->getCondType())));
	QList<QVariant> condData;
	condData << condType << "";

	const std::vector<GeoLib::Point*>* pnts = _project.getGEOObjects()->getPointVec(
			cond->getAssociatedGeometryName());
	if (pnts)
	{
		CondObjectListItem* cond_list = new CondObjectListItem(condData, parent, cond, pnts);
		parent->appendChild(cond_list);
		emit conditionAdded(this, parent->getItem()->getProcessType(), cond->getCondType());
		return cond_list;
	}

	return nullptr;
}

vtkPolyDataAlgorithm* ProcessModel::vtkSource(const FiniteElement::ProcessType pcs_type,
		const FEMCondition::CondType cond_type)
{
	ProcessItem* processParent = this->getProcessParent(pcs_type);
	if (processParent)
	{
		CondObjectListItem* condParent = this->getCondParent(processParent, cond_type);
		if (condParent)
			return condParent->vtkSource();
	}
	return nullptr;
}

void ProcessModel::replaceCondition(const QModelIndex &idx, FEMCondition* condition)
{
	// remove old condition
	this->getItem(idx)->parentItem()->removeChildren(this->getItem(idx)->row(), 1);
	//add new condition
	this->addCondition(condition);
}

void ProcessModel::updateModel()
{
	const std::vector<FEMCondition*> cnd_vec = _project.getConditions();
	for (auto it(cnd_vec.cbegin()); it != cnd_vec.cend(); ++it)
		// if Condition is not yet added to GUI, do it now
		if (!this->getCondition((*it)->getProcessType(), (*it)->getAssociatedGeometryName(), (*it)->getCondType(), (*it)->getGeomType(), (*it)->getGeoName())) 
			addConditionItem(*it);
}
