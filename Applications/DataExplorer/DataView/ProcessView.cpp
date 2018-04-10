/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <QFileDialog>
#include <QMenu>

#include "Applications/DataHolderLib/FemCondition.h"

#include "CondItem.h"
#include "ProcessModel.h"
#include "ProcessVarItem.h"
#include "ProcessView.h"
//#include "FEMConditionSetupDialog.h"
#include "SelectMeshDialog.h"

ProcessView::ProcessView(QWidget* parent) : QTreeView(parent) {}

void ProcessView::updateView()
{
    setAlternatingRowColors(true);
    resizeColumnToContents(0);
    setColumnWidth(1, 50);
    setColumnWidth(2, 50);
}

void ProcessView::on_Clicked(QModelIndex idx)
{
    qDebug("%d, %d", idx.parent().row(), idx.row());
}

void ProcessView::selectionChanged(const QItemSelection& selected,
                                   const QItemSelection& deselected)
{
    if (!selected.isEmpty())
    {
        emit clearConditionView();
        const QModelIndex idx = *(selected.indexes().begin());

        if (idx.parent().isValid())  // not root node
        {
            CondItem const*const item = dynamic_cast<CondItem*>(
                static_cast<ProcessModel*>(this->model())->getItem(idx));
            if (item != nullptr)
                emit conditionSelected(item->getCondition());
        }
        else
        {
            CondItem const*const item = dynamic_cast<CondItem*>(static_cast<ProcessModel*>(this->model())->getItem(idx)->child(0));
            if (item != nullptr)
                emit processVarSelected(item->getCondition());
        }
    }
    emit itemSelectionChanged(selected, deselected);
    QTreeView::selectionChanged(selected, deselected);
}

void ProcessView::contextMenuEvent(QContextMenuEvent* event)
{
    Q_UNUSED(event);

    const QModelIndex idx(this->selectionModel()->currentIndex());
    QMenu menu;

    if (isProcessVarItem(idx))
    {
        // QAction* saveCondAction  = menu.addAction("Save FEM Conditions...");
        QAction* removeProcessVarAction =
            menu.addAction("Remove process variable");
        // connect(saveCondAction, SIGNAL(triggered()), this,
        // SLOT(saveConditions()));
        connect(removeProcessVarAction, SIGNAL(triggered()),
                this, SLOT(removeProcessVar()));
    }
    else if (isConditionItem(idx))
    {
        QAction* removeCondAction = menu.addAction("Remove condition");
        connect(removeCondAction, SIGNAL(triggered()),
                this, SLOT(removeCondition()));
        // QAction* editCondAction = menu.addAction("Edit condition");
        // connect(editCondAction, SIGNAL(triggered()), this,
        // SLOT(editCondition()));  else  editCondAction->setEnabled(false);
    }

    menu.exec(event->globalPos());
}

void ProcessView::removeCondition()
{
    CondItem* item = dynamic_cast<CondItem*>(
        static_cast<ProcessModel*>(this->model())
            ->getItem(this->selectionModel()->currentIndex()));
    if (item)
    {
        emit clearConditionView();
        emit conditionRemoved(
            QString::fromStdString(item->getCondition()->getProcessVarName()),
            item->getName());
    }
}

void ProcessView::removeProcessVar()
{
    ProcessVarItem* item = dynamic_cast<ProcessVarItem*>(
        static_cast<ProcessModel*>(this->model())
            ->getItem(this->selectionModel()->currentIndex()));
    if (item)
    {
        emit clearConditionView();
        emit processVarRemoved(item->getName());
    }
}
/*
void ProcessView::editCondition()
{
    CondItem* item =
dynamic_cast<CondItem*>(static_cast<ProcessModel*>(this->model())->getItem(this->selectionModel()->currentIndex()));

    if (item)
    {
        FEMConditionSetupDialog dlg(*(item->getItem()));
        connect(&dlg, SIGNAL(createFEMCondition(std::vector<FEMCondition*>)),
this, SLOT(replaceCondition(std::vector<FEMCondition*>))); dlg.exec();
    }
}

void ProcessView::replaceCondition(std::vector<FEMCondition*> conditions)
{
    static_cast<ProcessModel*>(this->model())->replaceCondition(this->selectionModel()->currentIndex(),
conditions[0]); this->reset();
}

void ProcessView::saveConditions()
{
    emit saveConditionsRequested();
}
*/
bool ProcessView::isProcessVarItem(const QModelIndex& idx) const
{
    return (dynamic_cast<ProcessVarItem*>(
        static_cast<ProcessModel*>(this->model())->getItem(idx)) != nullptr);
}

bool ProcessView::isConditionItem(const QModelIndex& idx) const
{
    return (dynamic_cast<CondItem*>(
        static_cast<ProcessModel*>(this->model())->getItem(idx)) != nullptr);
}
