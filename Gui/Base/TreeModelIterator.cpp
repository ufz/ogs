/**
 * \file TreeModelIterator.cpp
 * 23/6/2010 LB Initial implementation
 *
 * Implementation of TreeModelIterator
 */

// ** INCLUDES **
#include "TreeModelIterator.h"

#include "TreeItem.h"
#include "TreeModel.h"

TreeModelIterator::TreeModelIterator( TreeModel* model )
	: _current(NULL), _model(model)
{
	if (_model->rootItem()->childCount() > 0)
	{
		_current = _model->rootItem();
		next(_current);
		//_parentIndex.push(0);
		//_currentIndex = 0;
	}
}

TreeItem* TreeModelIterator::operator*() const
{
	return _current;
}

TreeModelIterator& TreeModelIterator::operator++()
{
	if (_current)
		_current = next(_current);

	return *this;
}

TreeItem* TreeModelIterator::next( const TreeItem* current )
{
	if (!current)
		return NULL;

	TreeItem* next = NULL;
	if (current->childCount())
	{
		// walk the child
		_parentIndex.push(_currentIndex);
		_currentIndex = 0;
		next = current->child(0);
	}
	else
	{
		// walk the sibling
		TreeItem* parent = current->parentItem();
		next = parent ? parent->child(_currentIndex + 1)
		       : _model->rootItem()->child(_currentIndex + 1);
		while (!next && parent)
		{
			// if we had no sibling walk up the parent
			// and try the sibling of that
			parent = parent->parentItem();
			_currentIndex = _parentIndex.pop();
			next = parent ? parent->child(_currentIndex + 1)
			       : _model->rootItem()->child(_currentIndex + 1);
		}
		if (next)
			++(_currentIndex);
	}
	return next;
}
