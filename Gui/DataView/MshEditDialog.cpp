/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file MshEditDialog.cpp
 *
 * Created on 2010-11-09 by Karsten Rink
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

MshEditDialog::MshEditDialog(const MeshLib::Mesh* mesh, QDialog* parent)
	: QDialog(parent), _msh(mesh), _noDataDeleteBox(NULL), 
	  _nLayerLabel (new QLabel("Please specify the number of layers to add:")),  
	  _selectLabel(NULL),
	  _layerEdit (new QLineEdit("0")),
	  _nextButton (new QPushButton("Next"))
{
	setupUi(this);

	this->gridLayoutLayerMapping->addWidget(_nLayerLabel, 0, 0, 1, 2);
	this->gridLayoutLayerMapping->addWidget(_layerEdit, 0, 2);
	this->gridLayoutLayerMapping->addWidget(_nextButton, 0, 3);
	connect(_nextButton, SIGNAL(pressed()), this, SLOT(nextButtonPressed()));
}

MshEditDialog::~MshEditDialog()
{
	
	delete _nLayerLabel;
	delete _selectLabel;
	delete _layerEdit;
	delete _nextButton;
	delete _noDataDeleteBox;

	for (int i = 0; i < _labels.size(); i++)
	{
		delete _labels[i];
		delete _edits[i];
		delete _buttons[i];
	}
	
}

void MshEditDialog::nextButtonPressed()
{
	_layerEdit->setEnabled(false);
	_nextButton->setEnabled(false);
	const size_t nLayers = _layerEdit->text().toInt();
	const QString selectText = (nLayers>0) ?
		"Please specify a raster file for mapping each layer:" :
		"Please specify which rasterfile surface mapping:";
	_selectLabel = new QLabel(selectText);
	_selectLabel->setMargin(20);
	this->gridLayoutLayerMapping->addWidget(_selectLabel, 1, 0, 1, 4);

	this->gridLayoutLayerMapping->setMargin(10);
	this->gridLayoutLayerMapping->setColumnMinimumWidth(2,10);
	this->gridLayoutLayerMapping->setColumnStretch(0, 80);
	this->gridLayoutLayerMapping->setColumnStretch(1, 200);
	this->gridLayoutLayerMapping->setColumnStretch(2, 10);
	
	for (size_t i = 0; i <= nLayers+1; i++)
	{
		QString text("");
		if (i==0) text="Surface";
		else if (i>nLayers) text="Layer" + QString::number(nLayers) + "-Bottom";
		else text="Layer" + QString::number(i) + "-Top";
		QLabel* label = new QLabel(text);
		QLineEdit* edit = new QLineEdit();
		QPushButton* button = new QPushButton("...");

		_labels.push_back(label);
		_edits.push_back(edit);
		_buttons.push_back(button);
		_fileButtonMap.insert(button, edit);
		connect(button, SIGNAL(clicked()), this, SLOT(getFileName()));

		this->gridLayoutLayerMapping->addWidget(_labels[i],   i+2, 0);
		this->gridLayoutLayerMapping->addWidget(_edits[i],    i+2, 1, 1, 2);
		this->gridLayoutLayerMapping->addWidget(_buttons[i],  i+2, 3);

		if (nLayers==0) break; // don't add bottom layer if mesh contains only surface
	}

	_noDataDeleteBox = new QCheckBox("Remove mesh nodes at NoData values");
	_noDataDeleteBox->setChecked(false);
	_noDataDeleteBox->setEnabled(false);
	if (nLayers == 0)
	{
		_noDataDeleteBox->setEnabled(true);
		this->gridLayoutLayerMapping->addWidget(_noDataDeleteBox, 4, 1);
	}
}

void MshEditDialog::accept()
{
	if (_labels.size()>0)
	{
		bool all_paths_set (true);
		if ((_labels.size()==1) && (_edits[0]->text().length()==0))
			all_paths_set = false;
		else
		{
			for (int i=1; i<_labels.size(); i++)
				if (_edits[i]->text().length()==0)
					all_paths_set = false;
		}

		if (all_paths_set)
		{
			const size_t nLayers = _layerEdit->text().toInt();
			MeshLib::Mesh* new_mesh (NULL);

			if (nLayers==0)
			{
				new_mesh = new MeshLib::Mesh(*_msh);
				const std::string imgPath ( this->_edits[0]->text().toStdString() );
				if (!imgPath.empty())
					MshLayerMapper::LayerMapping(new_mesh, imgPath, nLayers, 0, _noDataDeleteBox->isChecked());
			}
			else
			{
				new_mesh = MshLayerMapper::CreateLayers(_msh, nLayers, 100);

				for (size_t i = 0; i <= nLayers; i++)
				{
					const std::string imgPath ( this->_edits[i]->text().toStdString() );
					if (!imgPath.empty())
					{
						int result = MshLayerMapper::LayerMapping(new_mesh, imgPath, nLayers, i, _noDataDeleteBox->isChecked());
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

			if (new_mesh)
			{
				new_mesh->setName("NewMesh");
				emit mshEditFinished(new_mesh);
			}
			else
				OGSError::box("Error creating mesh");

			this->done(QDialog::Accepted);
		}
		else
			OGSError::box("Please specifiy raster files for all layers.");
	}
	else
		OGSError::box("Please specifiy the number of\n layers and press \"Next\"");
}

void MshEditDialog::reject()
{
	this->done(QDialog::Rejected);
}

void MshEditDialog::getFileName()
{
	QPushButton* button = dynamic_cast<QPushButton*>(this->sender());
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString filename = QFileDialog::getOpenFileName(this, "Select raster file to open",
	                                                settings.value("lastOpenedFileDirectory").toString(),
	                                                "ASCII raster files (*.asc);;All files (* *.*)");
	_fileButtonMap[button]->setText(filename);
	QDir dir = QDir(filename);
	settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
}
