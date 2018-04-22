/**
 * \file
 * \author Karsten Rink
 * \date   2010-12-13
 * \brief  Definition of the ProcessView class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QContextMenuEvent>
#include <QTreeView>

#include "Applications/DataHolderLib/FemCondition.h"

class ConditionModel;

/**
 * \brief A view for FEM-Conditions (Initial- & Boundary Conditions / Source
 * Terms) with a number of additional information such as Process Type,
 * Distribution, etc. \sa ConditionModel, CondItem
 */
class ProcessView final : public QTreeView
{
    Q_OBJECT

public:
    /// Constructor
    ProcessView(QWidget* parent = nullptr);

    /// Update the view to visualise changes made to the underlying data
    void updateView();

protected slots:
    /// Instructions if the selection of items in the view has changed.
    void selectionChanged(const QItemSelection& selected,
                          const QItemSelection& deselected);

private:
    /// Actions to be taken after a right mouse click is performed in the
    /// station view.
    void contextMenuEvent(QContextMenuEvent* e);
    bool isProcessVarItem(const QModelIndex& idx) const;
    bool isConditionItem(const QModelIndex& idx) const;

private slots:
    void on_Clicked(QModelIndex idx);
    // void editCondition();
    void removeCondition();
    void removeProcessVar();
    // void replaceCondition(std::vector<FEMCondition*> conditions);
    // void saveConditions();

signals:
    void itemSelectionChanged(QItemSelection const& selected,
                              QItemSelection const& deselected);
    void conditionRemoved(QString const&, QString const&);
    void processVarRemoved(QString const&);
    // void saveConditionsRequested();
    void clearConditionView();
    void processVarSelected(DataHolderLib::FemCondition* cond);
    void conditionSelected(DataHolderLib::FemCondition* cond);
};
