/**
 * \file MshItem.cpp
 * 17/05/2010 KR Initial implementation
 */

#include "MshItem.h"
#include "GridAdapter.h"
#include "VtkMeshSource.h"



/**
 * Constructor.
 * \param data The data associated with each column
 * \param parent The parent item in the tree
 * \param grid The mesh associated with this item
 */
MshItem::MshItem(const QList<QVariant> &data, TreeItem *parent, GridAdapter* grid)
: TreeItem(data, parent)
{
	_meshSource = VtkMeshSource::New();
	_meshSource->SetGrid(grid);
}

MshItem::~MshItem()
{
	_meshSource->Delete();
}
