/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-09
 * \brief  Definition of the MeshLayerEditDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
    MeshLayerEditDialog(const MeshLib::Mesh* mesh, QDialog* parent = nullptr);
    ~MeshLayerEditDialog(void) override;

private:
    void createMeshToolSelection();
    MeshLib::Mesh* createPrismMesh();
    MeshLib::Mesh* createTetMesh();

    const MeshLib::Mesh* _msh;
    unsigned _n_layers;
    QMap<QPushButton*, QLineEdit*> _fileButtonMap;
    QVector<QLineEdit*> _edits;

    QLineEdit* _layerEdit;
    QLineEdit* _noDataReplacementEdit;
    QLineEdit* _minThicknessEdit;
    QPushButton* _nextButton;
    QGroupBox* _layerBox;
    QGroupBox* _radioButtonBox;
    QRadioButton* _ogsMeshButton;
    QGridLayout* _layerSelectionLayout;
    bool _use_rasters;

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
