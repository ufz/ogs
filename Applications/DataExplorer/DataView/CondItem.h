// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

#include "Applications/DataHolderLib/FemCondition.h"
#include "Base/TreeItem.h"

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

    /// Returns the FEM Condition associated with the item.
    DataHolderLib::FemCondition* getCondition() const { return _cond; }

    QString const getName() const { return data(0).toString(); }

private:
    DataHolderLib::FemCondition* _cond;
};
