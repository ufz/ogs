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

#include <memory>

#include "Applications/DataHolderLib/FemCondition.h"
#include "TreeItem.h"

#include <QList>
#include <QVariant>

/**
 * \brief A TreeItem containing a boundary condition or source term
 * \sa TreeItem
 */
class CondItem : public TreeItem
{
public:
    /// Constructor
    CondItem(const QList<QVariant>& data, TreeItem* parent,
             DataHolderLib::FemCondition* cond)
        : TreeItem(data, parent), _cond(cond)
    {
    }

    ~CondItem() = default;

    /// Returns the FEM Condition associated with the item.
    DataHolderLib::FemCondition* getCondition() const { return _cond; }

    QString const getName() const { return data(0).toString(); }

private:
    DataHolderLib::FemCondition* _cond;
};
