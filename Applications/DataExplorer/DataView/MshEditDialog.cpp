/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-09
 * \brief  Implementation of the MshEditDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MshEditDialog.h"

#include "OGSError.h"
#include "StringTools.h"
#include "Mesh.h"

#include <QCheckBox>
#include <QFileDialog>
#include <QPushButton>
#include <QSettings>
#include <QLineEdit>
#include <QGroupBox>
#include <QRadioButton>
#include <QVBoxLayout>

MshEditDialog::MshEditDialog(const MeshLib::Mesh* mesh, QDialog* parent)
    : QDialog(parent),
      _msh(mesh),
      _noDataDeleteBox(nullptr),
      _nLayerLabel(new QLabel("Please specify the number of layers to add:")),
      _nLayerExplanation(new QLabel("(select \"0\" for surface mapping)")),
      _layerEdit(new QLineEdit("0")),
      _nextButton(new QPushButton("Next")),
      _layerBox(nullptr),
      _radioButtonBox(nullptr),
      _layerSelectionLayout(new QGridLayout),
      _radiobuttonLayout(new QVBoxLayout),
      _selectButton1(new QRadioButton("Add layers based on raster files")),
      _selectButton2(new QRadioButton("Add layers with static thickness")),
      _n_layers(0),
      _use_rasters(true)
{
    setupUi(this);

    this->gridLayoutLayerMapping->addWidget(_nLayerLabel, 0, 0);
    this->gridLayoutLayerMapping->addWidget(_nLayerExplanation, 1, 0);
    this->gridLayoutLayerMapping->addWidget(_layerEdit, 0, 1);
    this->gridLayoutLayerMapping->addWidget(_nextButton, 0, 2);
    connect(_nextButton, SIGNAL(pressed()), this, SLOT(nextButtonPressed()));
}

MshEditDialog::~MshEditDialog()
{
    for (int i = 0; i < _labels.size(); ++i)
    {
        delete _labels[i];
        delete _edits[i];
    }
    for (int i = 0; i < _buttons.size(); ++i)
        delete _buttons[i];

    delete _nLayerLabel;
    delete _nLayerExplanation;
    delete _layerEdit;
    delete _nextButton;
    delete _noDataDeleteBox;
    delete _radiobuttonLayout;
    delete _layerSelectionLayout;
    delete _layerBox;
    delete _selectButton1;
    delete _selectButton2;
}

void MshEditDialog::nextButtonPressed()
{
    _layerEdit->setEnabled(false);
    _nextButton->setEnabled(false);
    _nLayerExplanation->setText("");
    _n_layers = _layerEdit->text().toInt();

    // configure group box + layout (will be needed in the next step)
    _layerBox = new QGroupBox;
    this->_layerSelectionLayout->setMargin(10);
    this->_layerSelectionLayout->setColumnMinimumWidth(2,10);
    this->_layerSelectionLayout->setColumnStretch(0, 80);
    this->_layerSelectionLayout->setColumnStretch(1, 200);
    this->_layerSelectionLayout->setColumnStretch(2, 10);

    if (_n_layers > 0)
    {
        _radioButtonBox = new QGroupBox;
        _radiobuttonLayout->addWidget(_selectButton1);
        _radiobuttonLayout->addWidget(_selectButton2);
        _radioButtonBox->setLayout(_radiobuttonLayout);
        gridLayoutLayerMapping->addWidget(_radioButtonBox, 2, 0, 1, 3);
        // add an empty line to better arrange the following information
        gridLayoutLayerMapping->addWidget(_nLayerExplanation, 3, 0);
        connect(_selectButton1, SIGNAL(pressed()), this, SLOT(createWithRasters()));
        connect(_selectButton2, SIGNAL(pressed()), this, SLOT(createStatic()));
    }
    else
        this->createWithRasters();


}

void MshEditDialog::createWithRasters()
{
    // _use_rasters=true is needed for this, this is the default setting however
    //this->_radioButtonBox->setEnabled(false);
    const QString selectText = (_n_layers>0) ?
            "Please specify a raster file for mapping each layer:" :
            "Please specify rasterfile for surface mapping:";
    this->_layerBox->setTitle(selectText);

    for (unsigned i = 0; i <= _n_layers+1; ++i)
    {
        QString text("");
        if (i==0) text = "Surface";
        else if (i>_n_layers) text = "Layer" + QString::number(_n_layers) + "-Bottom";
        else text="Layer" + QString::number(i) + "-Top";
        QLineEdit* edit (new QLineEdit());
        QPushButton* button (new QPushButton("..."));

        this->_labels.push_back(new QLabel(text));
        this->_edits.push_back(edit);
        this->_buttons.push_back(button);
        this->_fileButtonMap.insert(button, edit);
        connect(button, SIGNAL(clicked()), this, SLOT(getFileName()));

        this->_layerSelectionLayout->addWidget(_labels[i],  i, 0);
        this->_layerSelectionLayout->addWidget(_edits[i],   i, 1);
        this->_layerSelectionLayout->addWidget(_buttons[i], i, 2);

        // don't add bottom layer if mesh contains only surface
        if (this->_n_layers==0) break;
    }
    this->_layerBox->setLayout(this->_layerSelectionLayout);
    this->_noDataDeleteBox = new QCheckBox("Remove mesh nodes at NoData values");
    this->_noDataDeleteBox->setChecked(false);
    this->_noDataDeleteBox->setEnabled(false);
    if (this->_n_layers == 0)
    {
        this->_noDataDeleteBox->setEnabled(true);
        this->_layerSelectionLayout->addWidget(_noDataDeleteBox, 1, 1);
    }
    this->gridLayoutLayerMapping->addWidget(_layerBox, 4, 0, 1, 3);
}

void MshEditDialog::createStatic()
{
    this->_use_rasters = false;
    this->_radioButtonBox->setEnabled(false);
    this->_layerBox->setTitle("Please specify a thickness or each layer");

    for (unsigned i = 0; i < this->_n_layers; ++i)
    {
        QString text("Layer" + QString::number(i) + "-Thickness");
        _labels.push_back(new QLabel(text));
        _edits.push_back(new QLineEdit());
        this->_layerSelectionLayout->addWidget(_labels[i],  i, 0);
        this->_layerSelectionLayout->addWidget(_edits[i],   i, 1);
    }
    this->_layerBox->setLayout(this->_layerSelectionLayout);
    this->gridLayoutLayerMapping->addWidget(_layerBox, 5, 0, 1, 3);
}

void MshEditDialog::accept()
{
    if (this->_labels.size()>0)
    {
        bool all_paths_set (true);
        if ((_n_layers==0) && _use_rasters && (_edits[0]->text().length()==0))
            all_paths_set = false;
        else
        {
            int start_idx = (_use_rasters) ? 1:0;
            for (int i=start_idx; i<_labels.size(); ++i)
                if (_edits[i]->text().length()==0)
                    all_paths_set = false;
        }

        if (all_paths_set)
        {
            int result(1);
            const unsigned nLayers = _layerEdit->text().toInt();
            MeshLib::Mesh* new_mesh(nullptr);

            if (nLayers==0)
            {
                new_mesh = new MeshLib::Mesh(*_msh);
                const std::string imgPath ( this->_edits[0]->text().toStdString() );
                if (!imgPath.empty())
                    result = MshLayerMapper::LayerMapping(new_mesh, imgPath, nLayers, 0, _noDataDeleteBox->isChecked());
            }
            else
            {
                std::vector<float> layer_thickness(_n_layers);
                for (unsigned i=0; i<nLayers; ++i)
                    layer_thickness[i] = (_use_rasters) ? 100 : this->_edits[i]->text().toFloat();

                new_mesh = MshLayerMapper::CreateLayers(_msh, layer_thickness);

                if (_use_rasters)
                {
                    for (unsigned i=0; i<=nLayers; ++i)
                    {
                        const std::string imgPath ( this->_edits[i+1]->text().toStdString() );
                        if (!imgPath.empty())
                        {
                            result = MshLayerMapper::LayerMapping(new_mesh, imgPath, nLayers, i, _noDataDeleteBox->isChecked());
                            if (result==0) break;
                        }
                    }
                    if (this->_edits[0]->text().length()>0)
                    {
                        MeshLib::Mesh* final_mesh = MshLayerMapper::blendLayersWithSurface(new_mesh, nLayers, this->_edits[0]->text().toStdString());
                        delete new_mesh;
                        new_mesh = final_mesh;
                    }
                }
            }

            if (new_mesh)
                emit mshEditFinished(new_mesh);

            if (!new_mesh || result==0)
                OGSError::box("Error creating mesh");

            this->done(QDialog::Accepted);
        }
        else
            OGSError::box("Please specifiy raster files for all layers.");
    }
    else
        OGSError::box("Please specifiy the number and\n type of layers and press \"Next\"");
}

void MshEditDialog::reject()
{
    this->done(QDialog::Rejected);
}

void MshEditDialog::getFileName()
{
    QPushButton* button = dynamic_cast<QPushButton*>(this->sender());
    QSettings settings;
    QString filename = QFileDialog::getOpenFileName(this, "Select raster file to open",
                                                    settings.value("lastOpenedFileDirectory").toString(),
                                                    "ASCII raster files (*.asc);;All files (* *.*)");
    _fileButtonMap[button]->setText(filename);
    QDir dir = QDir(filename);
    settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
}
