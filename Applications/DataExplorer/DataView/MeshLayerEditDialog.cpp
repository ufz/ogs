/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-09
 * \brief  Implementation of the MeshLayerEditDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshGenerators/LayeredVolume.h"
#include "MeshLayerEditDialog.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "OGSError.h"
#include "StringTools.h"
#include "Mesh.h"

#include "TetGenInterface.h"

#include <QCheckBox>
#include <QFileInfo>
#include <QFileDialog>
#include <QPushButton>
#include <QSettings>
#include <QLineEdit>
#include <QGroupBox>
#include <QRadioButton>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QTime>

MeshLayerEditDialog::MeshLayerEditDialog(const MeshLib::Mesh* mesh, QDialog* parent)
	: QDialog(parent), _msh(mesh), _n_layers(0),
	  _nLayerExplanation (new QLabel(this)),
	  _layerEdit (new QLineEdit("0", this)),
	  _noDataReplacementEdit(nullptr),
	  _nextButton (new QPushButton("Next", this)),
	  _layerBox (nullptr),
	  _radioButtonBox (nullptr),
	  _ogsMeshButton (nullptr),
	  _layerSelectionLayout (new QGridLayout(_layerBox)),
	  _use_rasters(true)
{
	setupUi(this);

	_nLayerExplanation->setText("(select \"0\" for surface mapping)");
	this->gridLayoutLayerMapping->addWidget(new QLabel("Please specify the number of layers to add:", this), 0, 0);
	this->gridLayoutLayerMapping->addWidget(_nLayerExplanation, 1, 0);
	this->gridLayoutLayerMapping->addWidget(_layerEdit, 0, 1);
	this->gridLayoutLayerMapping->addWidget(_nextButton, 0, 2);
	_layerEdit->setValidator(new QIntValidator(0,999,_layerEdit));
	connect(_nextButton, SIGNAL(pressed()), this, SLOT(nextButtonPressed()));

	// configure group box + layout
	this->_layerSelectionLayout->setMargin(10);
	this->_layerSelectionLayout->setColumnMinimumWidth(2,10);
	this->_layerSelectionLayout->setColumnStretch(0, 80);
	this->_layerSelectionLayout->setColumnStretch(1, 200);
	this->_layerSelectionLayout->setColumnStretch(2, 10);
}

MeshLayerEditDialog::~MeshLayerEditDialog()
{
}

void MeshLayerEditDialog::nextButtonPressed()
{
	_layerEdit->setEnabled(false);
	_nextButton->setEnabled(false);
	_nLayerExplanation->setText("");
	_n_layers  = static_cast<unsigned>(_layerEdit->text().toInt());

	if (_n_layers == 0)
		this->createWithRasters();
	else
	{
		QVBoxLayout* _radiobuttonLayout (new QVBoxLayout(_radioButtonBox));
		QRadioButton* selectButton1 (new QRadioButton("Add layers based on raster files", _radioButtonBox));
		QRadioButton* selectButton2 (new QRadioButton("Add layers with static thickness", _radioButtonBox));
		_radioButtonBox = new QGroupBox(this);
		_radiobuttonLayout->addWidget(selectButton1);
		_radiobuttonLayout->addWidget(selectButton2);
		_radioButtonBox->setLayout(_radiobuttonLayout);
		gridLayoutLayerMapping->addWidget(_radioButtonBox, 2, 0, 1, 3);
		// add an empty line to better arrange the following information
		gridLayoutLayerMapping->addWidget(_nLayerExplanation, 3, 0);
		connect(selectButton1, SIGNAL(pressed()), this, SLOT(createWithRasters()));
		connect(selectButton2, SIGNAL(pressed()), this, SLOT(createStatic()));
	}
}

void MeshLayerEditDialog::createWithRasters()
{
	// _use_rasters=true is needed for this, this is the default setting however
	if (_n_layers>0)
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
		QLineEdit* edit (new QLineEdit(this));
		QPushButton* button (new QPushButton("...", _layerBox));

		this->_edits.push_back(edit);
		this->_fileButtonMap.insert(button, edit);
		connect(button, SIGNAL(clicked()), this, SLOT(getFileName()));

		this->_layerSelectionLayout->addWidget(new QLabel(text, _layerBox),  i, 0);
		this->_layerSelectionLayout->addWidget(_edits[i],   i, 1);
		this->_layerSelectionLayout->addWidget(button, i, 2);

		// don't add bottom layer if mesh contains only surface
		if (this->_n_layers==0) break;
	}
	if (this->_n_layers == 0)
	{
		QLabel* noDataReplacementLabel = new QLabel("Set NoData values to ", _layerBox);
		_noDataReplacementEdit = new QLineEdit("0.0", _layerBox);
		_noDataReplacementEdit->setValidator(new QDoubleValidator(_noDataReplacementEdit));

		this->_layerSelectionLayout->addWidget(noDataReplacementLabel, 1, 0);
		this->_layerSelectionLayout->addWidget(_noDataReplacementEdit, 1, 1);
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
	QGroupBox* meshToolSelectionBox (new QGroupBox(this));
	meshToolSelectionBox->setTitle("Select output element type");
	QHBoxLayout* meshToolSelectionLayout (new QHBoxLayout(meshToolSelectionBox));
	_ogsMeshButton = new QRadioButton("Prisms", meshToolSelectionBox);
	QRadioButton* tetgenMeshButton = new QRadioButton("Tetrahedra", meshToolSelectionBox);
	meshToolSelectionLayout->addWidget(_ogsMeshButton);
	meshToolSelectionLayout->addWidget(tetgenMeshButton);
	meshToolSelectionBox->setLayout(meshToolSelectionLayout);
	_ogsMeshButton->setChecked(true);
	gridLayoutLayerMapping->addWidget(meshToolSelectionBox, 5, 0, 1, 3);
}

MeshLib::Mesh* MeshLayerEditDialog::createPrismMesh()
{
    const unsigned nLayers = _layerEdit->text().toInt();
    std::vector<float> layer_thickness;
    for (unsigned i=0; i<nLayers; ++i)
    {
        // "100" is just a default size to have any value for extruding 2D elements.
        // The actual mapping based on a raster file will be performed later.
        const float thickness = (_use_rasters) ? 100 : (this->_edits[i]->text().toFloat());
        layer_thickness.push_back(thickness);
    }

    MeshLayerMapper mapper;
    MeshLib::Mesh* new_mesh (nullptr);

    QTime myTimer0;
    myTimer0.start();
    if (_use_rasters)
    {
        std::vector<std::string> raster_paths;
        for (int i=nLayers; i>=0; --i)
            raster_paths.push_back(this->_edits[i]->text().toStdString());
        if (mapper.createLayers(*_msh, raster_paths))
            new_mesh= mapper.getMesh("SubsurfaceMesh");
    }
    else
        new_mesh = mapper.createStaticLayers(*_msh, layer_thickness);
    INFO("Mesh construction time: %d ms.", myTimer0.elapsed());

    return new_mesh;
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
	MeshLib::Mesh* tg_mesh (nullptr);
	QTime myTimer0;
	myTimer0.start();

	if (_use_rasters)
	{
		std::vector<std::string> raster_paths;
		for (int i=nLayers; i>=0; --i)
			raster_paths.push_back(this->_edits[i]->text().toStdString());
		LayeredVolume lv;
		if (lv.createLayers(*_msh, raster_paths))
			tg_mesh = lv.getMesh("SubsurfaceMesh");

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
		MeshLayerMapper const mapper;
		tg_mesh = mapper.createStaticLayers(*_msh, layer_thickness);
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

	const unsigned nLayers = _layerEdit->text().toInt();
	MeshLib::Mesh* new_mesh (NULL);

	if (nLayers==0)
	{
		new_mesh = new MeshLib::Mesh(*_msh);
		const std::string imgPath ( this->_edits[0]->text().toStdString() );
		const double noDataReplacementValue = this->_noDataReplacementEdit->text().toDouble();
		MeshLayerMapper const mapper;
		if (!mapper.layerMapping(*new_mesh, imgPath, noDataReplacementValue))
		{
			delete new_mesh;
			return;
		}
	}
	else
	{
		if (_ogsMeshButton->isChecked())
			new_mesh = this->createPrismMesh();
		else
			new_mesh = this->createTetMesh();
	}

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
	QPushButton* button = dynamic_cast<QPushButton*>(this->sender());
	QSettings settings;
	QString filename = QFileDialog::getOpenFileName(this, "Select raster file to open",
	                                                settings.value("lastOpenedRasterFileDirectory").toString(),
	                                                "ASCII raster files (*.asc);;All files (* *.*)");
	_fileButtonMap[button]->setText(filename);
	QFileInfo fi(filename);
	settings.setValue("lastOpenedRasterFileDirectory", fi.absolutePath());
}
