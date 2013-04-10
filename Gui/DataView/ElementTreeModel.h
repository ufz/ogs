/**
 * \file
 * \author Karsten Rink
 * \date   2011-05-10
 * \brief  Definition of the ElementTreeModel class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENTTREEMODEL_H
#define ELEMENTTREEMODEL_H

#include "TreeModel.h"

class vtkUnstructuredGridAlgorithm;

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

	vtkUnstructuredGridAlgorithm const*const getSource() const { return _mesh_source; };

public slots:
	/// Clears the tree.
	void clearView();

	/// Extracts information of the element with the given index from the given grid.
	void setElement(vtkUnstructuredGridAlgorithm const*const grid, const std::size_t elem_index);

private:
	vtkUnstructuredGridAlgorithm const* _mesh_source;
};

#endif //ELEMENTTREEMODEL_H
