/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-09
 * \brief  Implementation of the MeshLayerEditDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshGenerators/LayeredVolume.h"
#include "MeshLayerEditDialog.h"

#include <logog/include/logog.hpp>

#include "OGSError.h"
#include "StringTools.h"
#include "Mesh.h"

#include "Applications/FileIO/AsciiRasterInterface.h"
#include "TetGenInterface.h"

#include <QCheckBox>
#include <QFileInfo>
#include <QFileDialog>
#include <QPushButton>
#include <QSettings>
#include <QLineEdit>
#include <QGroupBox>
#include <QRadioButton>
#include <QGridLayout>
#include <QVBoxLayout>
#include <QTime>

MeshLayerEditDialog::MeshLayerEditDialog(const MeshLib::Mesh* mesh, QDialog* parent)
    : QDialog(parent), _msh(mesh), _n_layers(0),
      _layerEdit (new QLineEdit("1", this)),
      _noDataReplacementEdit(nullptr),
      _minThicknessEdit(nullptr),
      _nextButton (new QPushButton("Next", this)),
      _layerBox (nullptr),
      _radioButtonBox (nullptr),
      _ogsMeshButton (nullptr),
      _layerSelectionLayout (new QGridLayout(_layerBox)),
      _use_rasters(true)
{
    setupUi(this);

    this->gridLayoutLayerMapping->addWidget(new QLabel("Please specify the number of layers to add:", this), 0, 0);
    this->gridLayoutLayerMapping->addWidget(_layerEdit, 0, 1);
    this->gridLayoutLayerMapping->addWidget(_nextButton, 0, 2);
    _layerEdit->setValidator(new QIntValidator(1,999,_layerEdit));
    connect(_nextButton, SIGNAL(pressed()), this, SLOT(nextButtonPressed()));

    // configure group box + layout
    this->_layerSelectionLayout->setMargin(10);
    this->_layerSelectionLayout->setColumnMinimumWidth(2,10);
    this->_layerSelectionLayout->setColumnStretch(0, 80);
    this->_layerSelectionLayout->setColumnStretch(1, 200);
    this->_layerSelectionLayout->setColumnStretch(2, 10);
}

MeshLayerEditDialog::~MeshLayerEditDialog() = default;

void MeshLayerEditDialog::nextButtonPressed()
{
    _n_layers  = static_cast<unsigned>(_layerEdit->text().toInt());

    if (_n_layers < 1)
    {
        OGSError::box("Add the number of layers to add (at least 1)");
        return;
    }

    _layerEdit->setEnabled(false);
    _nextButton->setEnabled(false);

    auto* _radiobuttonLayout(new QVBoxLayout(_radioButtonBox));
    QRadioButton* selectButton1 (new QRadioButton("Add layers based on raster files", _radioButtonBox));
    QRadioButton* selectButton2 (new QRadioButton("Add layers with static thickness", _radioButtonBox));
    _radioButtonBox = new QGroupBox(this);
    _radiobuttonLayout->addWidget(selectButton1);
    _radiobuttonLayout->addWidget(selectButton2);
    _radioButtonBox->setLayout(_radiobuttonLayout);
    gridLayoutLayerMapping->addWidget(_radioButtonBox, 2, 0, 1, 3);
    connect(selectButton1, SIGNAL(pressed()), this, SLOT(createWithRasters()));
    connect(selectButton2, SIGNAL(pressed()), this, SLOT(createStatic()));
}

void MeshLayerEditDialog::createWithRasters()
{
    this->_use_rasters = true;
    this->_radioButtonBox->setEnabled(false);
    const QString selectText = (_n_layers>0) ?
            "Please specify a raster file for mapping each layer:" :
            "Please specify raster file for surface mapping:";
    this->_layerBox = new QGroupBox(this);
    this->_layerBox->setTitle(selectText);

    for (unsigned i = 0; i <= _n_layers; ++i)
    {
        QString text("");
        if (i==0) text = "Surface";
        else if (i == _n_layers) text = "Layer" + QString::number(_n_layers) + "-Bottom";
        else text="Layer" + QString::number(i+1) + "-Top";
        auto* edit(new QLineEdit(this));
        QPushButton* button (new QPushButton("...", _layerBox));

        this->_edits.push_back(edit);
        this->_fileButtonMap.insert(button, edit);
        connect(button, SIGNAL(clicked()), this, SLOT(getFileName()));

        this->_layerSelectionLayout->addWidget(new QLabel(text, _layerBox),  i, 0);
        this->_layerSelectionLayout->addWidget(_edits[i],   i, 1);
        this->_layerSelectionLayout->addWidget(button, i, 2);
    }
    this->_layerBox->setLayout(this->_layerSelectionLayout);
    this->gridLayoutLayerMapping->addWidget(_layerBox, 4, 0, 1, 3);
    if (this->_n_layers > 0)
        this->createMeshToolSelection();
}

void MeshLayerEditDialog::createStatic()
{
    this->_use_rasters = false;
    this->_radioButtonBox->setEnabled(false);
    this->_layerBox = new QGroupBox(this);
    this->_layerBox->setTitle("Please specify a thickness or each layer");

    for (unsigned i = 0; i < this->_n_layers; ++i)
    {
        QString text("Layer" + QString::number(i) + "-Thickness");
        QLineEdit* staticLayerEdit = new QLineEdit("10", this);
        staticLayerEdit->setValidator(new QDoubleValidator(staticLayerEdit));
        _edits.push_back(staticLayerEdit);
        this->_layerSelectionLayout->addWidget(new QLabel(text, _layerBox),  i, 0);
        this->_layerSelectionLayout->addWidget(_edits[i],   i, 1);
    }
    this->_layerBox->setLayout(this->_layerSelectionLayout);
    this->gridLayoutLayerMapping->addWidget(_layerBox, 4, 0, 1, 3);
    this->createMeshToolSelection();
}

void MeshLayerEditDialog::createMeshToolSelection()
{
    auto* meshToolSelectionBox(new QGroupBox(this));
    meshToolSelectionBox->setTitle("Output element type");
    auto* meshToolSelectionLayout(new QGridLayout(meshToolSelectionBox));
    _ogsMeshButton = new QRadioButton("Prisms", meshToolSelectionBox);
    QRadioButton* tetgenMeshButton = new QRadioButton("Tetrahedra", meshToolSelectionBox);
    tetgenMeshButton->setFixedWidth(150);
    QLabel* minThicknessLabel = new QLabel(meshToolSelectionBox);
    minThicknessLabel->setText("Minimum thickness of layers:");
    _minThicknessEdit = new QLineEdit(meshToolSelectionBox);
    _minThicknessEdit->setText("1.0");
    auto* min_thickness_validator =
        new QDoubleValidator(0, 1000000, 15, _minThicknessEdit);
    _minThicknessEdit->setValidator(min_thickness_validator);
    _minThicknessEdit->setMaximumWidth(100);
    _minThicknessEdit->setFixedWidth(100);
    meshToolSelectionLayout->addWidget(_ogsMeshButton, 0, 0);
    meshToolSelectionLayout->addWidget(tetgenMeshButton, 0, 1);
    meshToolSelectionLayout->addWidget(minThicknessLabel, 1, 0);
    meshToolSelectionLayout->addWidget(_minThicknessEdit, 1, 1);
    meshToolSelectionBox->setLayout(meshToolSelectionLayout);
    _ogsMeshButton->setChecked(true);

    gridLayoutLayerMapping->addWidget(meshToolSelectionBox, 5, 0, 1, 3);
}

MeshLib::Mesh* MeshLayerEditDialog::createPrismMesh()
{
    const unsigned nLayers = _layerEdit->text().toInt();

    MeshLib::MeshLayerMapper mapper;

    QTime myTimer0;
    myTimer0.start();
    if (_use_rasters)
    {
        float minimum_thickness (_minThicknessEdit->text().toFloat());
        if (minimum_thickness <= 0) minimum_thickness = std::numeric_limits<float>::epsilon();
        std::vector<std::string> raster_paths;
        for (int i=nLayers; i>=0; --i)
            raster_paths.push_back(this->_edits[i]->text().toStdString());

        auto const rasters = FileIO::readRasters(raster_paths);
        if (rasters && mapper.createLayers(*_msh, *rasters, minimum_thickness))
        {
            INFO("Mesh construction time: %d ms.", myTimer0.elapsed());
            return mapper.getMesh("SubsurfaceMesh").release();
        }
        return nullptr;
    }

    std::vector<float> layer_thickness;
    for (unsigned i = 0; i < nLayers; ++i)
        layer_thickness.push_back(this->_edits[i]->text().toFloat());
    INFO("Mesh construction time: %d ms.", myTimer0.elapsed());
    return mapper.createStaticLayers(*_msh, layer_thickness);
}

MeshLib::Mesh* MeshLayerEditDialog::createTetMesh()
{
    QSettings settings;
    QString filename = QFileDialog::getSaveFileName(this, "Write TetGen input file to",
                                                    settings.value("lastOpenedTetgenFileDirectory").toString(),
                                                    "TetGen Geometry (*.smesh)");
    if (filename.isEmpty())
        return nullptr;

    const unsigned nLayers = _layerEdit->text().toInt();
    MeshLib::Mesh* tg_mesh(nullptr);
    QTime myTimer0;
    myTimer0.start();

    if (_use_rasters)
    {
        float minimum_thickness (_minThicknessEdit->text().toFloat());
        if (minimum_thickness <= 0) minimum_thickness = std::numeric_limits<float>::epsilon();
        std::vector<std::string> raster_paths;
        for (int i=nLayers; i>=0; --i)
            raster_paths.push_back(this->_edits[i]->text().toStdString());
        LayeredVolume lv;

        auto const rasters = FileIO::readRasters(raster_paths);
        if (rasters && lv.createLayers(*_msh, *rasters, minimum_thickness))
            tg_mesh = lv.getMesh("SubsurfaceMesh").release();

        if (tg_mesh)
        {
            std::vector<MeshLib::Node> tg_attr (lv.getAttributePoints());
            FileIO::TetGenInterface tetgen_interface;
            tetgen_interface.writeTetGenSmesh(filename.toStdString(), *tg_mesh, tg_attr);
        }
    }
    else
    {
        std::vector<float> layer_thickness;
        for (unsigned i=0; i<nLayers; ++i)
            layer_thickness.push_back(this->_edits[i]->text().toFloat());
        tg_mesh = MeshLib::MeshLayerMapper::createStaticLayers(*_msh, layer_thickness);
        std::vector<MeshLib::Node> tg_attr;
        FileIO::TetGenInterface tetgen_interface;
        tetgen_interface.writeTetGenSmesh(filename.toStdString(), *tg_mesh, tg_attr);
    }
    INFO("Mesh construction time: %d ms.", myTimer0.elapsed());

    return tg_mesh;
}

void MeshLayerEditDialog::accept()
{
    if (this->_edits.isEmpty())
    {
        OGSError::box("Please specifiy the number and\n type of layers and press \"Next\"");
        return;
    }

    bool all_paths_set (true);
    if (_n_layers==0)
    {
        if (_edits[0]->text().isEmpty())
            all_paths_set = false;
    }
    else
    {
        int start_idx = (_use_rasters) ? 1:0;
        for (int i=start_idx; i<_edits.size(); ++i)
            if (_edits[i]->text().isEmpty())
                all_paths_set = false;
    }

    if (!all_paths_set)
    {
        OGSError::box("Please specifiy raster files for all layers.");
        return;
    }

    MeshLib::Mesh* new_mesh (nullptr);
    if (_ogsMeshButton->isChecked())
        new_mesh = createPrismMesh();
    else
        new_mesh = createTetMesh();

    if (new_mesh)
        emit mshEditFinished(new_mesh);
    else
        OGSError::box("Error creating mesh");

    this->done(QDialog::Accepted);
}

void MeshLayerEditDialog::reject()
{
    this->done(QDialog::Rejected);
}

void MeshLayerEditDialog::getFileName()
{
    auto* button = dynamic_cast<QPushButton*>(this->sender());
    QSettings settings;
    QString filename = QFileDialog::getOpenFileName(
        this, "Select raster file to open",
        settings.value("lastOpenedRasterFileDirectory").toString(),
        "ASCII raster files (*.asc);;All files (* *.*)");
    _fileButtonMap[button]->setText(filename);
    QFileInfo fi(filename);
    settings.setValue("lastOpenedRasterFileDirectory", fi.absolutePath());
}
