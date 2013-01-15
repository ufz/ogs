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

namespace MeshLib {
	class Mesh;
}

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
	void setElement(const MeshLib::Mesh* grid, const std::size_t elem_index);

private:
};

#endif //ELEMENTTREEMODEL_H
