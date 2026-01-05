// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base/TreeItem.h"

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
