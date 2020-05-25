/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-09
 * \brief  Definition of the MeshLayerEditDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ui_MeshLayerEdit.h"
#include <QDialog>
#include <QLineEdit>

#include "MeshGenerators/MeshLayerMapper.h"

class QPushButton;
class QCheckBox;
class QGroupBox;
class QGridLayout;
class QVBoxLayout;
class QRadioButton;

namespace MeshLib
{
class Mesh;
}

/**
 * \brief A dialog window for editing meshes in various ways
 */
class MeshLayerEditDialog : public QDialog, private Ui_MeshLayerEdit
{
    Q_OBJECT

public:
    explicit MeshLayerEditDialog(const MeshLib::Mesh* mesh,
                                 QDialog* parent = nullptr);
    ~MeshLayerEditDialog() override;

private:
    void createMeshToolSelection();
    MeshLib::Mesh* createPrismMesh();
    MeshLib::Mesh* createTetMesh();

    const MeshLib::Mesh* msh_;
    unsigned n_layers_;
    QMap<QPushButton*, QLineEdit*> fileButtonMap_;
    QVector<QLineEdit*> edits_;

    QLineEdit* layerEdit_;
    QLineEdit* noDataReplacementEdit_;
    QLineEdit* minThicknessEdit_;
    QPushButton* nextButton_;
    QGroupBox* layerBox_;
    QGroupBox* radioButtonBox_;
    QRadioButton* ogsMeshButton_;
    QGridLayout* layerSelectionLayout_;
    bool use_rasters_;

private slots:
    void getFileName();

    void nextButtonPressed();

    void createWithRasters();

    void createStatic();

    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override;

signals:
    void mshEditFinished(MeshLib::Mesh*);
};
