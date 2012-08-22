/**
 * \file TreeModelIterator.h
 * 23/6/2010 LB Initial implementation
 *
 */

#ifndef TREEMODELITERATOR_H
#define TREEMODELITERATOR_H

// ** INCLUDES **
#include <QStack>

class TreeModel;
class TreeItem;

/**
 * \brief TreeModelIterator provides a way to iterate over TreeItems in a TreeModel.
 * Usage: \code
 * TreeModelIterator it(model);
 * while(*it)
 * {
 *		QVariant var = (*it)->data(0);
 *		std::cout << var.toString().toStdString() << std::endl;
 *		++it;
 * } \endcode
 */
class TreeModelIterator
{
public:
	/// \brief Constructor. Provide a tree model to iterate over.
	TreeModelIterator(TreeModel* model);

	/// \brief Dereferencing the iterator to retrieve the current TreeItem.
	/// Returns NULL if the iterator is at the end.
	TreeItem* operator* () const;

	/// \brief Advance the iterator to the next TreeItem.
	TreeModelIterator& operator++ ();

private:
	/// \brief The current TreeItem.
	TreeItem* _current;

	/// \brief The current child index.
	int _currentIndex;

	/// \brief Stack to save the child indices of the parent TreeItems.
	QStack<int> _parentIndex;

	/// \brief The model to iterate over.
	TreeModel* _model;

	/// \brief The traversal implementation.
	TreeItem* next(const TreeItem* current);
};

#endif // TREEMODELITERATOR_H
