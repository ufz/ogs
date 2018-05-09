/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-16
 * \brief  Definition of the MeshElementRemovalDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ui_MeshElementRemoval.h"

#include <QDialog>
#include <array>

#include "MeshLib/PropertyVector.h"


namespace DataHolderLib
{
class Project;
}

namespace MeshLib {
class Mesh;
}

/**
 * \brief A dialog window for settung up a database connection
 */
class MeshElementRemovalDialog : public QDialog, private Ui_MeshElementRemoval
{
    Q_OBJECT

public:
    MeshElementRemovalDialog(DataHolderLib::Project const& project,
                             QDialog* parent = nullptr);
    ~MeshElementRemovalDialog(void) override;

private slots:
    void on_boundingBoxCheckBox_toggled(bool is_checked);
    void on_elementTypeCheckBox_toggled(bool is_checked);
    void on_scalarArrayCheckBox_toggled(bool is_checked);
    void on_insideButton_toggled(bool is_checked);
    void on_meshNameComboBox_currentIndexChanged(int idx);
    void on_scalarArrayComboBox_currentIndexChanged(int idx);
    void on_xMinEdit_textChanged() { aabb_edits[0] = true; }
    void on_xMaxEdit_textChanged() { aabb_edits[1] = true; }
    void on_yMinEdit_textChanged() { aabb_edits[2] = true; }
    void on_yMaxEdit_textChanged() { aabb_edits[3] = true; }
    void on_zMinEdit_textChanged() { aabb_edits[4] = true; }
    void on_zMaxEdit_textChanged() { aabb_edits[5] = true; }
    void accept() override;
    void reject() override;

private:
    std::size_t addScalarArrays(MeshLib::Mesh const& mesh) const;
    void enableScalarArrayWidgets(bool enable) const;

    template <typename T>
    void setRangeValues(MeshLib::PropertyVector<T> const& vec);

    void toggleScalarEdits(bool outside) const;

    DataHolderLib::Project const& _project;
    unsigned _currentIndex, _aabbIndex, _scalarIndex;
    std::array<bool, 6> aabb_edits;

signals:
    void meshAdded(MeshLib::Mesh* mesh);
};
