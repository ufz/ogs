/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Applications/DataHolderLib/FemCondition.h"
#include "TreeModel.h"

/**
 * \brief A model for the display of information from boundary conditions and
 * source terms. \sa TreeModel, FemConditionView, TreeItem
 */
class FemConditionModel : public TreeModel
{
    Q_OBJECT

public:
    FemConditionModel(QObject* parent = nullptr);
    ~FemConditionModel() override;

    int columnCount() const { return 2; }

public slots:
    /// Clears the tree.
    void clearView();

    /// Displays information on a boundary condition or source term
    void setFemCondition(DataHolderLib::FemCondition* cond);

    /// Displays information on a process variable
    void setProcessVariable(DataHolderLib::FemCondition* cond);

private:
};
