/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-09
 * \brief  Implementation of the MeshLayerEditDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshGenerators/LayeredVolume.h"
#include "MeshLayerEditDialog.h"

#include "BaseLib/Logging.h"

#include "OGSError.h"
#include "StringTools.h"
#include "Mesh.h"

#include "Applications/FileIO/AsciiRasterInterface.h"
#include "TetGenInterface.h"

#include <QCheckBox>
#include <QDoubleValidator>
#include <QFileDialog>
#include <QFileInfo>
#include <QGridLayout>
#include <QGroupBox>
#include <QIntValidator>
#include <QLineEdit>
#include <QPushButton>
#include <QRadioButton>
#include <QSettings>
#include <QElapsedTimer>
#include <QVBoxLayout>

MeshLayerEditDialog::MeshLayerEditDialog(const MeshLib::Mesh* mesh, QDialog* parent)
    : QDialog(parent), msh_(mesh), n_layers_(0),
      layerEdit_ (new QLineEdit("1", this)),
      noDataReplacementEdit_(nullptr),
      minThicknessEdit_(nullptr),
      nextButton_ (new QPushButton("Next", this)),
      layerBox_ (nullptr),
      radioButtonBox_ (nullptr),
      ogsMeshButton_ (nullptr),
      layerSelectionLayout_ (new QGridLayout(layerBox_)),
      use_rasters_(true)
{
    setupUi(this);

    this->gridLayoutLayerMapping->addWidget(new QLabel("Please specify the number of layers to add:", this), 0, 0);
    this->gridLayoutLayerMapping->addWidget(layerEdit_, 0, 1);
    this->gridLayoutLayerMapping->addWidget(nextButton_, 0, 2);
    layerEdit_->setValidator(new QIntValidator(1,999,layerEdit_));
    connect(nextButton_, SIGNAL(pressed()), this, SLOT(nextButtonPressed()));

    // configure group box + layout
    this->layerSelectionLayout_->setMargin(10);
    this->layerSelectionLayout_->setColumnMinimumWidth(2,10);
    this->layerSelectionLayout_->setColumnStretch(0, 80);
    this->layerSelectionLayout_->setColumnStretch(1, 200);
    this->layerSelectionLayout_->setColumnStretch(2, 10);
}

MeshLayerEditDialog::~MeshLayerEditDialog() = default;

void MeshLayerEditDialog::nextButtonPressed()
{
    n_layers_  = static_cast<unsigned>(layerEdit_->text().toInt());

    if (n_layers_ < 1)
    {
        OGSError::box("Add the number of layers to add (at least 1)");
        return;
    }

    layerEdit_->setEnabled(false);
    nextButton_->setEnabled(false);

    auto* radiobuttonLayout_(new QVBoxLayout(radioButtonBox_));
    QRadioButton* selectButton1 (new QRadioButton("Add layers based on raster files", radioButtonBox_));
    QRadioButton* selectButton2 (new QRadioButton("Add layers with static thickness", radioButtonBox_));
    radioButtonBox_ = new QGroupBox(this);
    radiobuttonLayout_->addWidget(selectButton1);
    radiobuttonLayout_->addWidget(selectButton2);
    radioButtonBox_->setLayout(radiobuttonLayout_);
    gridLayoutLayerMapping->addWidget(radioButtonBox_, 2, 0, 1, 3);
    connect(selectButton1, SIGNAL(pressed()), this, SLOT(createWithRasters()));
    connect(selectButton2, SIGNAL(pressed()), this, SLOT(createStatic()));
}

void MeshLayerEditDialog::createWithRasters()
{
    this->use_rasters_ = true;
    this->radioButtonBox_->setEnabled(false);
    const QString selectText = (n_layers_>0) ?
            "Please specify a raster file for mapping each layer:" :
            "Please specify raster file for surface mapping:";
    this->layerBox_ = new QGroupBox(this);
    this->layerBox_->setTitle(selectText);

    for (unsigned i = 0; i <= n_layers_; ++i)
    {
        QString text("");
        if (i == 0)
        {
            text = "Surface";
        }
        else if (i == n_layers_)
        {
            text = "Layer" + QString::number(n_layers_) + "-Bottom";
        }
        else
        {
            text = "Layer" + QString::number(i + 1) + "-Top";
        }
        auto* edit(new QLineEdit(this));
        QPushButton* button (new QPushButton("...", layerBox_));

        this->edits_.push_back(edit);
        this->fileButtonMap_.insert(button, edit);
        connect(button, SIGNAL(clicked()), this, SLOT(getFileName()));

        this->layerSelectionLayout_->addWidget(new QLabel(text, layerBox_),  i, 0);
        this->layerSelectionLayout_->addWidget(edits_[i],   i, 1);
        this->layerSelectionLayout_->addWidget(button, i, 2);
    }
    this->layerBox_->setLayout(this->layerSelectionLayout_);
    this->gridLayoutLayerMapping->addWidget(layerBox_, 4, 0, 1, 3);
    if (this->n_layers_ > 0)
    {
        this->createMeshToolSelection();
    }
}

void MeshLayerEditDialog::createStatic()
{
    this->use_rasters_ = false;
    this->radioButtonBox_->setEnabled(false);
    this->layerBox_ = new QGroupBox(this);
    this->layerBox_->setTitle("Please specify a thickness or each layer");

    for (unsigned i = 0; i < this->n_layers_; ++i)
    {
        QString text("Layer" + QString::number(i) + "-Thickness");
        QLineEdit* staticLayerEdit = new QLineEdit("10", this);
        staticLayerEdit->setValidator(new QDoubleValidator(staticLayerEdit));
        edits_.push_back(staticLayerEdit);
        this->layerSelectionLayout_->addWidget(new QLabel(text, layerBox_),  i, 0);
        this->layerSelectionLayout_->addWidget(edits_[i],   i, 1);
    }
    this->layerBox_->setLayout(this->layerSelectionLayout_);
    this->gridLayoutLayerMapping->addWidget(layerBox_, 4, 0, 1, 3);
    this->createMeshToolSelection();
}

void MeshLayerEditDialog::createMeshToolSelection()
{
    auto* meshToolSelectionBox(new QGroupBox(this));
    meshToolSelectionBox->setTitle("Output element type");
    auto* meshToolSelectionLayout(new QGridLayout(meshToolSelectionBox));
    ogsMeshButton_ = new QRadioButton("Prisms", meshToolSelectionBox);
    QRadioButton* tetgenMeshButton = new QRadioButton("Tetrahedra", meshToolSelectionBox);
    tetgenMeshButton->setFixedWidth(150);
    auto* minThicknessLabel = new QLabel(meshToolSelectionBox);
    minThicknessLabel->setText("Minimum thickness of layers:");
    minThicknessEdit_ = new QLineEdit(meshToolSelectionBox);
    minThicknessEdit_->setText("1.0");
    auto* min_thickness_validator =
        new QDoubleValidator(0, 1000000, 15, minThicknessEdit_);
    minThicknessEdit_->setValidator(min_thickness_validator);
    minThicknessEdit_->setMaximumWidth(100);
    minThicknessEdit_->setFixedWidth(100);
    meshToolSelectionLayout->addWidget(ogsMeshButton_, 0, 0);
    meshToolSelectionLayout->addWidget(tetgenMeshButton, 0, 1);
    meshToolSelectionLayout->addWidget(minThicknessLabel, 1, 0);
    meshToolSelectionLayout->addWidget(minThicknessEdit_, 1, 1);
    meshToolSelectionBox->setLayout(meshToolSelectionLayout);
    ogsMeshButton_->setChecked(true);

    gridLayoutLayerMapping->addWidget(meshToolSelectionBox, 5, 0, 1, 3);
}

MeshLib::Mesh* MeshLayerEditDialog::createPrismMesh()
{
    const unsigned nLayers = layerEdit_->text().toInt();

    MeshLib::MeshLayerMapper mapper;

    QElapsedTimer myTimer0;
    myTimer0.start();
    if (use_rasters_)
    {
        float minimum_thickness (minThicknessEdit_->text().toFloat());
        if (minimum_thickness <= 0)
        {
            minimum_thickness = std::numeric_limits<float>::epsilon();
        }
        std::vector<std::string> raster_paths;
        for (int i = nLayers; i >= 0; --i)
        {
            raster_paths.push_back(this->edits_[i]->text().toStdString());
        }

        auto const rasters = FileIO::readRasters(raster_paths);
        if (rasters && mapper.createLayers(*msh_, *rasters, minimum_thickness))
        {
            INFO("Mesh construction time: {:d} ms.", myTimer0.elapsed());
            return mapper.getMesh("SubsurfaceMesh").release();
        }
        return nullptr;
    }

    std::vector<float> layer_thickness;
    for (unsigned i = 0; i < nLayers; ++i)
    {
        layer_thickness.push_back(this->edits_[i]->text().toFloat());
    }
    INFO("Mesh construction time: {:d} ms.", myTimer0.elapsed());
    return MeshLib::MeshLayerMapper::createStaticLayers(*msh_, layer_thickness);
}

MeshLib::Mesh* MeshLayerEditDialog::createTetMesh()
{
    QSettings settings;
    QString filename = QFileDialog::getSaveFileName(this, "Write TetGen input file to",
                                                    settings.value("lastOpenedTetgenFileDirectory").toString(),
                                                    "TetGen Geometry (*.smesh)");
    if (filename.isEmpty())
    {
        return nullptr;
    }

    const unsigned nLayers = layerEdit_->text().toInt();
    MeshLib::Mesh* tg_mesh(nullptr);
    QElapsedTimer myTimer0;
    myTimer0.start();

    if (use_rasters_)
    {
        float minimum_thickness (minThicknessEdit_->text().toFloat());
        if (minimum_thickness <= 0)
        {
            minimum_thickness = std::numeric_limits<float>::epsilon();
        }
        std::vector<std::string> raster_paths;
        for (int i = nLayers; i >= 0; --i)
        {
            raster_paths.push_back(this->edits_[i]->text().toStdString());
        }
        LayeredVolume lv;

        auto const rasters = FileIO::readRasters(raster_paths);
        if (rasters && lv.createLayers(*msh_, *rasters, minimum_thickness))
        {
            tg_mesh = lv.getMesh("SubsurfaceMesh").release();
        }

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
        for (unsigned i = 0; i < nLayers; ++i)
        {
            layer_thickness.push_back(this->edits_[i]->text().toFloat());
        }
        tg_mesh = MeshLib::MeshLayerMapper::createStaticLayers(*msh_, layer_thickness);
        std::vector<MeshLib::Node> tg_attr;
        FileIO::TetGenInterface tetgen_interface;
        tetgen_interface.writeTetGenSmesh(filename.toStdString(), *tg_mesh, tg_attr);
    }
    INFO("Mesh construction time: {:d} ms.", myTimer0.elapsed());

    return tg_mesh;
}

void MeshLayerEditDialog::accept()
{
    if (this->edits_.isEmpty())
    {
        OGSError::box(
            "Please specifiy the number and\n type of layers and press 'Next'");
        return;
    }

    bool all_paths_set (true);
    if (n_layers_==0)
    {
        if (edits_[0]->text().isEmpty())
        {
            all_paths_set = false;
        }
    }
    else
    {
        int start_idx = (use_rasters_) ? 1:0;
        for (int i = start_idx; i < edits_.size(); ++i)
        {
            if (edits_[i]->text().isEmpty())
            {
                all_paths_set = false;
            }
        }
    }

    if (!all_paths_set)
    {
        OGSError::box("Please specifiy raster files for all layers.");
        return;
    }

    MeshLib::Mesh* new_mesh (nullptr);
    if (ogsMeshButton_->isChecked())
    {
        new_mesh = createPrismMesh();
    }
    else
    {
        new_mesh = createTetMesh();
    }

    if (new_mesh)
    {
        emit mshEditFinished(new_mesh);
    }
    else
    {
        OGSError::box("Error creating mesh");
    }

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
    fileButtonMap_[button]->setText(filename);
    QFileInfo fi(filename);
    settings.setValue("lastOpenedRasterFileDirectory", fi.absolutePath());
}
