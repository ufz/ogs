/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSITEM_H
#define PROCESSITEM_H

#include "TreeItem.h"

/**
 * \brief A TreeItem representing process variable information.
 * \sa TreeItem
 */
class ProcessVarItem : public TreeItem
{
public:
    /// Constructor
    ProcessVarItem(const QList<QVariant>& data, TreeItem* parent)
        : TreeItem(data, parent)
    {
    }

    ~ProcessVarItem() {}

    QString const& getName() const { return data(0).toString(); }
};

#endif  // PROCESSITEM_H
