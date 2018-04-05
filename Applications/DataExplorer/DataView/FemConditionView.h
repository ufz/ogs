/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <QTreeView>

#include "Applications/DataHolderLib/FemCondition.h"

// class vtkUnstructuredGridAlgorithm;

/**
 *    A TreeView to display information of FEM conditions.
 */
class FemConditionView : public QTreeView
{
    Q_OBJECT

public:
    /// Constructor
    FemConditionView(QWidget* parent = nullptr);

public slots:
    void updateView();

protected slots:
    /// Is called when the selection of this view changes.
    void selectionChanged(const QItemSelection& selected,
                          const QItemSelection& deselected) override;

signals:
};
