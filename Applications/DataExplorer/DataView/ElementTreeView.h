/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-09
 * \brief  Definition of the ElementTreeView class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <QTreeView>

class vtkUnstructuredGridAlgorithm;

/**
 *    A TreeView to display mesh element properties.
 */
class ElementTreeView : public QTreeView
{
    Q_OBJECT

public:
    /// Constructor
    ElementTreeView(QWidget* parent = nullptr);

public slots:
    void updateView();

protected slots:
    /// Is called when the selection of this view changes.
    void selectionChanged(const QItemSelection& selected,
                          const QItemSelection& deselected) override;

signals:
    void nodeSelected(vtkUnstructuredGridAlgorithm const*const, unsigned, bool);
    void removeSelectedMeshComponent();
};
