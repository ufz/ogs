/**
 * \file
 * \author Karsten Rink
 * \date   2010-10-18
 * \brief  Implementation of the ProcessModel class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessModel.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// ** INCLUDES **
#include "GeoLib/GEOObjects.h"
#include "GeoLib/GeoObject.h"
#include "Applications/DataHolderLib/FemCondition.h"

#include "ProcessVarItem.h"
#include "CondItem.h"
#include "GeoType.h"

#include <QFileInfo>
#include <vtkPolyDataAlgorithm.h>

ProcessModel::ProcessModel(DataHolderLib::Project &project, QObject* parent)
: TreeModel(parent), _project(project)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Name" << "Value" << "" << "" << "";
	_rootItem = new TreeItem(rootData, nullptr);
}

ProcessModel::~ProcessModel()
{
}

int ProcessModel::columnCount(QModelIndex const& parent) const
{
	Q_UNUSED(parent)

	return 2;
}

void ProcessModel::addConditionItem(DataHolderLib::FemCondition* cond, ProcessVarItem* parent)
{
    QList<QVariant> item_data;

    item_data << QString::fromStdString(cond->getParamName())
                << QString::fromStdString(cond->getConditionClassStr());

    CondItem* cond_item = new CondItem(item_data, parent, cond);
    parent->appendChild(cond_item);
}

void ProcessModel::addCondition(DataHolderLib::FemCondition* condition)
{
    QString name(QString::fromStdString(condition->getProcessVarName()));
    ProcessVarItem* process_var (getProcessVarItem(name));
    if (process_var == nullptr)
        process_var = addProcessVar(name);
        addConditionItem(condition, process_var);
}

void ProcessModel::addBoundaryConditions(std::vector<std::unique_ptr<DataHolderLib::BoundaryCondition>> const& conditions)
{
    for (size_t i = 0; i < conditions.size(); i++)
        addCondition(conditions[i].get());
}

void ProcessModel::addSourceTerms(std::vector<std::unique_ptr<DataHolderLib::SourceTerm>> const& conditions)
{
    for (size_t i = 0; i < conditions.size(); i++)
        addCondition(conditions[i].get());
}

ProcessVarItem* ProcessModel::addProcessVar(QString const& name)
{
    beginResetModel();
    QList<QVariant> process_var_data;
    process_var_data << QVariant(name) << "";
    ProcessVarItem* process_var = new ProcessVarItem(process_var_data, _rootItem);
    _rootItem->appendChild(process_var);
    endResetModel();
    return process_var;
}

ProcessVarItem* ProcessModel::getProcessVarItem(QString const& name) const
{
    int n_children (_rootItem->childCount());
    for (int i=0; i<n_children; ++i)
    {
        ProcessVarItem* item (dynamic_cast<ProcessVarItem*>(_rootItem->child(i)));
        if (item != nullptr && item->getName() == name)
            return item;
    }
    return nullptr;
}

void ProcessModel::removeCondition(ProcessVarItem* process_var, QString const& param_name)
{
    int n_conditions = process_var->childCount();
    for (int i = 0; i<n_conditions; ++i)
    {
        CondItem const*const cond = dynamic_cast<CondItem*>(process_var->child(i));
        if (cond->getCondition()->getParamName() != param_name.toStdString())
            continue;

        process_var->removeChildren(i, 1);
        return;
    }
}

void ProcessModel::removeCondition(QString const& process_var, QString const& param)
{
    beginResetModel();
    ProcessVarItem* pv_item(getProcessVarItem(process_var));
    if (pv_item == nullptr)
        return;

    removeCondition(pv_item, param);
    _project.removeBoundaryCondition(process_var.toStdString(), param.toStdString());
    _project.removeSourceTerm(process_var.toStdString(), param.toStdString());
    endResetModel();
}

void ProcessModel::removeProcessVariable(QString const& name)
{
    beginResetModel();
    ProcessVarItem* pv_item (getProcessVarItem(name));
    if (pv_item == nullptr)
        return;

    int n_conds = pv_item->childCount();
    for (int i=n_conds-1; i>=0; --i)
        removeCondition(pv_item, static_cast<CondItem*>(pv_item->child(i))->getName());

    _project.removePrimaryVariable(name.toStdString());
    int idx = pv_item->row();
    _rootItem->removeChildren(idx, 1);
    endResetModel();
}

void ProcessModel::clearModel()
{
    int n_process_vars = _rootItem->childCount();
    for (int i = n_process_vars; i>=0; --i)
    {
        ProcessVarItem* pv_item = dynamic_cast<ProcessVarItem*>(_rootItem->child(i));
        removeProcessVariable(pv_item->getName());
    }
}

void ProcessModel::updateModel()
{
    addBoundaryConditions(_project.getBoundaryConditions());
    addSourceTerms(_project.getSourceTerms());
}

