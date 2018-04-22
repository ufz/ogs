/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "TreeItem.h"

/**
 * \brief A TreeItem representing process variable information.
 * \sa TreeItem
 */
class ProcessVarItem final : public TreeItem
{
public:
    /// Constructor
    ProcessVarItem(const QList<QVariant>& data, TreeItem* parent)
        : TreeItem(data, parent)
    {
    }

    QString getName() const { return data(0).toString(); }
};
