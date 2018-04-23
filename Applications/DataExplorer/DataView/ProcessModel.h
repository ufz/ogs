/**
 * \file
 * \author Karsten Rink
 * \date   2010-10-18
 * \brief  Definition of the ProcessModel class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Applications/DataHolderLib/Project.h"
#include "TreeModel.h"

class FEMCondition;
class ProcessVarItem;
class CondObjectListItem;
class vtkPolyDataAlgorithm;

namespace GeoLib
{
class GeoObject;
}

/**
 * \brief A model implementing a tree structure for process-relevant information
 * such as process types, FEM-Conditions (BCs, ICs, STs), etc. as a
 * double-linked list. \sa TreeModel, ProcessView, TreeItem, CondObjectListItem
 */
class ProcessModel final : public TreeModel
{
    Q_OBJECT

public:
    ProcessModel(DataHolderLib::Project& project, QObject* parent = nullptr);

    int columnCount(const QModelIndex& parent = QModelIndex()) const;

    /// Returns the vtk source object for the specified subtree of a process
    /// with the given name.
    // vtkPolyDataAlgorithm* vtkSource(const FiniteElement::ProcessType
    // pcs_type, const FEMCondition::CondType cond_type);

public slots:
    /// Adds vector of Boundary Conditions to the model.
    void addBoundaryConditions(
        std::vector<std::unique_ptr<DataHolderLib::BoundaryCondition>> const&
            conditions);

    /// Adds vector of Source Terms to the model.
    void addSourceTerms(
        std::vector<std::unique_ptr<DataHolderLib::SourceTerm>> const&
            conditions);

    /// Adds a single FEM Conditions to the model
    void addCondition(DataHolderLib::FemCondition* condition);

    /// Adds a process to the model
    ProcessVarItem* addProcessVar(QString const& name);

    /// Removes FEMConditions from the the model
    void removeCondition(QString const& process_var, QString const& param);

    /// Removes a process variable incl all associated conditions from the model
    void removeProcessVariable(QString const& name);

    /// Removes the complete content from the model
    void clearModel();

    /// Updates the model based on the ProjectData-object
    void updateModel();

private:
    /// Adds a new FEM condition to the condition tree model.
    void addConditionItem(DataHolderLib::FemCondition* cond,
                          ProcessVarItem* parent);

    ProcessVarItem* getProcessVarItem(QString const& process_var_name) const;

    /// Removes FEMConditions from the the model
    void removeCondition(ProcessVarItem* parent, QString const& param_name);

    DataHolderLib::Project& _project;

signals:
};
