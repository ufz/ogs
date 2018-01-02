/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-09
 * \brief  Implementation of the ElementTreeView class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include "ElementTreeView.h"

#include "ElementTreeModel.h"
#include "TreeItem.h"

#include <QModelIndex>

ElementTreeView::ElementTreeView( QWidget* parent)
: QTreeView(parent)
{
}

void ElementTreeView::updateView()
{
    setAlternatingRowColors(true);
    setColumnWidth(0,125);
    std::size_t nColumns =
        (this->model() != nullptr) ? this->model()->columnCount() : 0;
    for (std::size_t i = 1; i < nColumns; i++)
        resizeColumnToContents(i);
    this->expandAll();
}

void ElementTreeView::selectionChanged( const QItemSelection &selected, const QItemSelection &deselected )
{
    Q_UNUSED(deselected);
    if (!selected.isEmpty())
    {
        emit removeSelectedMeshComponent();
        const QModelIndex idx = *(selected.indexes().begin());

        if (idx.parent().isValid()) // not root node
            if (idx.parent().parent().isValid()) // not property node
            {
                const TreeItem* tree_item = static_cast<TreeModel*>(this->model())->getItem(idx);
                const unsigned node_index = tree_item->data(0).toString().mid(5).toUInt();
                emit nodeSelected(static_cast<ElementTreeModel*>(this->model())->getSource(), node_index, false);
            }
    }
}
