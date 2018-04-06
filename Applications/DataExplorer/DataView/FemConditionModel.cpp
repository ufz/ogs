/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FemConditionModel.h"

#include "Applications/DataHolderLib/BoundaryCondition.h"
#include "Applications/DataHolderLib/SourceTerm.h"

#include "TreeItem.h"

/**
 * Constructor.
 */
FemConditionModel::FemConditionModel(QObject* parent) : TreeModel(parent)
{
    QList<QVariant> root_data;
    delete _rootItem;
    root_data << "Parameter"
              << "Value";
    _rootItem = new TreeItem(root_data, nullptr);
}

FemConditionModel::~FemConditionModel() = default;

void FemConditionModel::setFemCondition(DataHolderLib::FemCondition* cond)
{
    beginResetModel();
    this->clearView();

    QList<QVariant> cond_data;
    cond_data << QString::fromStdString(cond->getConditionClassStr()) + ":"
              << QString::fromStdString(cond->getParamName());
    TreeItem* cond_item = new TreeItem(cond_data, _rootItem);
    _rootItem->appendChild(cond_item);

    QList<QVariant> type_data;
    std::string type_str;
    if (cond->getConditionClassStr() == "Boundary Condition")
        type_str = DataHolderLib::BoundaryCondition::convertTypeToString(
            static_cast<DataHolderLib::BoundaryCondition*>(cond)->getType());
    else if (cond->getConditionClassStr() == "Source Term")
        type_str = DataHolderLib::SourceTerm::convertTypeToString(
            static_cast<DataHolderLib::SourceTerm*>(cond)->getType());
    type_data << "Type:" << QString::fromStdString(type_str);
    TreeItem* type_item = new TreeItem(type_data, cond_item);
    cond_item->appendChild(type_item);

    QString const obj_class_str =
        (cond->getBaseObjType() == DataHolderLib::BaseObjType::MESH)
            ? "Mesh:"
            : "Geometry:";
    QList<QVariant> obj_class_data;
    obj_class_data << "Set on " + obj_class_str
                   << QString::fromStdString(cond->getBaseObjName());
    TreeItem* obj_class_item = new TreeItem(obj_class_data, cond_item);
    cond_item->appendChild(obj_class_item);

    QString const obj_str =
        (cond->getBaseObjType() == DataHolderLib::BaseObjType::MESH)
            ? "Mesh array:"
            : "Geo-Object:";
    QString const name_str =
        (cond->getBaseObjType() == DataHolderLib::BaseObjType::MESH)
            ? QString::fromStdString(cond->getParamName())
            : QString::fromStdString(cond->getObjName());
    QList<QVariant> obj_data;
    obj_data << obj_str << name_str;
    TreeItem* obj_item = new TreeItem(obj_data, cond_item);
    cond_item->appendChild(obj_item);

    endResetModel();
}

void FemConditionModel::setProcessVariable(DataHolderLib::FemCondition* cond)
{
    beginResetModel();
    this->clearView();

    DataHolderLib::ProcessVariable const& var (cond->getProcessVar());

    QList<QVariant> pvar_data;
    pvar_data << "Process variable:" << QString::fromStdString(var.name);
    TreeItem* pvar_item = new TreeItem(pvar_data, _rootItem);
    _rootItem->appendChild(pvar_item);

    QList<QVariant> order_data;
    order_data << "Order:" << QString::number(var.order);
    TreeItem* order_item = new TreeItem(order_data, pvar_item);
    pvar_item->appendChild(order_item);

    QList<QVariant> comp_data;
    comp_data << "Number of components:" << QString::number(var.components);
    TreeItem* comp_item = new TreeItem(comp_data, pvar_item);
    pvar_item->appendChild(comp_item);

    endResetModel();
}

void FemConditionModel::clearView()
{
    beginResetModel();
    _rootItem->removeChildren(0, _rootItem->childCount());
    endResetModel();
}

