// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "FemConditionView.h"

#include <QModelIndex>

#include "CondItem.h"
#include "FemConditionModel.h"

FemConditionView::FemConditionView(QWidget* parent) : QTreeView(parent) {}

void FemConditionView::updateView()
{
    setAlternatingRowColors(true);
    setColumnWidth(0, 200);
    std::size_t nColumns =
        (this->model() != nullptr) ? this->model()->columnCount() : 0;
    for (std::size_t i = 1; i < nColumns; i++)
    {
        resizeColumnToContents(i);
    }
    this->expandAll();
}

void FemConditionView::selectionChanged(const QItemSelection& selected,
                                        const QItemSelection& deselected)
{
    Q_UNUSED(selected);
    Q_UNUSED(deselected);
}
