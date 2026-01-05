// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    explicit ElementTreeView(QWidget* parent = nullptr);

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
