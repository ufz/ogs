/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include "FemConditionView.h"

#include "CondItem.h"
#include "FemConditionModel.h"

#include <QModelIndex>

FemConditionView::FemConditionView(QWidget* parent) : QTreeView(parent) {}

void FemConditionView::updateView()
{
    setAlternatingRowColors(true);
    setColumnWidth(0, 200);
    std::size_t nColumns =
        (this->model() != nullptr) ? this->model()->columnCount() : 0;
    for (std::size_t i = 1; i < nColumns; i++)
        resizeColumnToContents(i);
    this->expandAll();
}

void FemConditionView::selectionChanged(const QItemSelection& selected,
                                        const QItemSelection& deselected)
{
    Q_UNUSED(deselected);
}
