/**
 * \file ElementTreeModel.h
 * 2011/05/10 KR Initial implementation
 */

#ifndef ELEMENTTREEMODEL_H
#define ELEMENTTREEMODEL_H

#include "GridAdapter.h"
#include "TreeModel.h"

/**
 * \brief A model for the display of information concerning element information implemented as a TreeModel.
 * \sa TreeModel, ElementTreeView, TreeItem
 */
class ElementTreeModel : public TreeModel
{
	Q_OBJECT

public:
	ElementTreeModel( QObject* parent = 0 );
	~ElementTreeModel();

public slots:
	/// Clears the tree.
	void clearView();

	/// Extracts information of the element with the given index from the given grid.
	void setElement(const GridAdapter* grid, const size_t elem_index);

private:
};

#endif //ELEMENTTREEMODEL_H
