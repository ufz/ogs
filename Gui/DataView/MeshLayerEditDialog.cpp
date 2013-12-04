/**
 * \file
 * \author Karsten Rink
 * \date   2010-11-09
 * \brief  Implementation of the MeshLayerEditDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#include "MeshLayerEditDialog.h"

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

MeshLayerEditDialog::MeshLayerEditDialog(const MeshLib::Mesh* mesh, QDialog* parent)
	: QDialog(parent), _msh(mesh), _n_layers(0),
	  _nLayerExplanation (new QLabel("(select \"0\" for surface mapping)")),
	  _layerEdit (new QLineEdit("0")),
	  _noDataReplacementEdit(nullptr),
	  _nextButton (new QPushButton("Next")),
	  _layerBox (nullptr),
	  _radioButtonBox (nullptr),
	  _layerSelectionLayout (new QGridLayout),
	  _use_rasters(true)
{
	setupUi(this);

	this->gridLayoutLayerMapping->addWidget(new QLabel("Please specify the number of layers to add:", this), 0, 0);
	this->gridLayoutLayerMapping->addWidget(_nLayerExplanation, 1, 0);
	this->gridLayoutLayerMapping->addWidget(_layerEdit, 0, 1);
	this->gridLayoutLayerMapping->addWidget(_nextButton, 0, 2);
	_layerEdit->setValidator(new QIntValidator(0,999,_layerEdit));
	connect(_nextButton, SIGNAL(pressed()), this, SLOT(nextButtonPressed()));

	// configure group box + layout
	_layerBox = new QGroupBox;
	this->_layerSelectionLayout->setMargin(10);
	this->_layerSelectionLayout->setColumnMinimumWidth(2,10);
	this->_layerSelectionLayout->setColumnStretch(0, 80);
	this->_layerSelectionLayout->setColumnStretch(1, 200);
	this->_layerSelectionLayout->setColumnStretch(2, 10);
}

MeshLayerEditDialog::~MeshLayerEditDialog()
{
	for (int i = 0; i < _edits.size(); ++i)
		delete _edits[i];

	delete _nLayerExplanation;
	delete _layerEdit;
	delete _noDataReplacementEdit;
	delete _nextButton;
	delete _radioButtonBox;
	delete _layerSelectionLayout;
	delete _layerBox;
}

void MeshLayerEditDialog::nextButtonPressed()
{
	_layerEdit->setEnabled(false);
	_nextButton->setEnabled(false);
	_nLayerExplanation->setText("");
	_n_layers  = static_cast<unsigned>(_layerEdit->text().toInt());

	if (_n_layers > 0)
	{
		QVBoxLayout* _radiobuttonLayout (new QVBoxLayout(_radioButtonBox));
		QRadioButton* selectButton1 (new QRadioButton("Add layers based on raster files", _radioButtonBox));
		QRadioButton* selectButton2 (new QRadioButton("Add layers with static thickness", _radioButtonBox));
		_radioButtonBox = new QGroupBox;
		_radiobuttonLayout->addWidget(selectButton1);
		_radiobuttonLayout->addWidget(selectButton2);
		_radioButtonBox->setLayout(_radiobuttonLayout);
		gridLayoutLayerMapping->addWidget(_radioButtonBox, 2, 0, 1, 3);
		// add an empty line to better arrange the following information
		gridLayoutLayerMapping->addWidget(_nLayerExplanation, 3, 0);
		connect(selectButton1, SIGNAL(pressed()), this, SLOT(createWithRasters()));
		connect(selectButton2, SIGNAL(pressed()), this, SLOT(createStatic()));
	}
	else
		this->createWithRasters();


}

void MeshLayerEditDialog::createWithRasters()
{
	// _use_rasters=true is needed for this, this is the default setting however
	if (_n_layers>0)
		this->_radioButtonBox->setEnabled(false);
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
}

void MeshLayerEditDialog::createStatic()
{
	this->_use_rasters = false;
	this->_radioButtonBox->setEnabled(false);
	this->_layerBox->setTitle("Please specify a thickness or each layer");

	for (unsigned i = 0; i < this->_n_layers; ++i)
	{
		QString text("Layer" + QString::number(i) + "-Thickness");
		QLineEdit* staticLayerEdit = new QLineEdit("10");
		staticLayerEdit->setValidator(new QDoubleValidator(staticLayerEdit));
		_edits.push_back(staticLayerEdit);
		this->_layerSelectionLayout->addWidget(new QLabel(text, _layerBox),  i, 0);
		this->_layerSelectionLayout->addWidget(_edits[i],   i, 1);
	}
	this->_layerBox->setLayout(this->_layerSelectionLayout);
	this->gridLayoutLayerMapping->addWidget(_layerBox, 5, 0, 1, 3);
}

void MeshLayerEditDialog::accept()
{
	if (this->_edits.size()>0)
	{
		bool all_paths_set (true);
		if ((_n_layers==0) && _use_rasters && (_edits[0]->text().length()==0))
			all_paths_set = false;
		else
		{
			int start_idx = (_use_rasters) ? 1:0;
			for (int i=start_idx; i<_edits.size(); ++i)
				if (_edits[i]->text().length()==0)
					all_paths_set = false;
		}

		if (all_paths_set)
		{
			int result(1);
			const unsigned nLayers = _layerEdit->text().toInt();
			MeshLib::Mesh* new_mesh (NULL);

			if (nLayers==0)
			{
				new_mesh = new MeshLib::Mesh(*_msh);
				const std::string imgPath ( this->_edits[0]->text().toStdString() );
				const double noDataReplacementValue = strtod(this->_noDataReplacementEdit->text().toStdString().c_str(),0);
				if (!imgPath.empty())
					result = MshLayerMapper::LayerMapping(new_mesh, imgPath, nLayers, 0, noDataReplacementValue);
			}
			else
			{
				std::vector<float> layer_thickness;
				for (unsigned i=0; i<nLayers; ++i)
				{
					float thickness = (_use_rasters) ? 100 : (this->_edits[i]->text().toFloat());
					if (thickness > std::numeric_limits<float>::epsilon())
						layer_thickness.push_back(thickness);
				}

				new_mesh = MshLayerMapper::CreateLayers(_msh, layer_thickness);

				if (_use_rasters)
				{
					for (unsigned i=0; i<=nLayers; ++i)
					{
						const std::string imgPath ( this->_edits[i+1]->text().toStdString() );
						const double noDataReplacement = (i==0) ? 0.0 : -9999.0;
						if (!imgPath.empty())
						{
							result = MshLayerMapper::LayerMapping(new_mesh, imgPath, nLayers, i, noDataReplacement);
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

void MeshLayerEditDialog::reject()
{
	this->done(QDialog::Rejected);
}

void MeshLayerEditDialog::getFileName()
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
