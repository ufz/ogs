// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Applications/DataHolderLib/FemCondition.h"
#include "Applications/DataExplorer/Base/TreeModel.h"

/**
 * \brief A model for the display of information from boundary conditions and
 * source terms. \sa TreeModel, FemConditionView, TreeItem
 */
class FemConditionModel : public TreeModel
{
    Q_OBJECT

public:
    explicit FemConditionModel(QObject* parent = nullptr);

    int columnCount(
        const QModelIndex& /*parent*/ = QModelIndex()) const override
    {
        return 2;
    }

public slots:
    /// Clears the tree.
    void clearView();

    /// Displays information on a boundary condition or source term
    void setFemCondition(DataHolderLib::FemCondition* cond);

    /// Displays information on a process variable
    void setProcessVariable(DataHolderLib::FemCondition* cond);
};
