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

MshEditDialog::MshEditDialog(const MeshLib::Mesh* mesh, QDialog* parent)
	: QDialog(parent), _msh(mesh), _noDataDeleteBox(NULL)
{
	setupUi(this);
/* TODO6
	this->gridLayoutLayerMapping->setMargin(5);
	this->gridLayoutLayerMapping->setColumnMinimumWidth(2,10);
	this->gridLayoutLayerMapping->setColumnStretch(0, 80);
	this->gridLayoutLayerMapping->setColumnStretch(1, 200);
	this->gridLayoutLayerMapping->setColumnStretch(2, 10);

	size_t nLayers = mesh->getNumberOfMeshLayers();

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

		this->gridLayoutLayerMapping->addWidget(_labels[i],   i, 0);
		this->gridLayoutLayerMapping->addWidget(_edits[i],    i, 1);
		this->gridLayoutLayerMapping->addWidget(_buttons[i],  i, 2);

		if (nLayers==0) break; // don't add bottom layer if mesh contains only surface
	}

	_noDataDeleteBox = new QCheckBox("Remove mesh nodes at NoData values");
	_noDataDeleteBox->setChecked(false);
	_noDataDeleteBox->setEnabled(false);
	if (nLayers == 0)
	{
		_noDataDeleteBox->setEnabled(true);
		this->gridLayoutLayerMapping->addWidget(_noDataDeleteBox, 2, 1);
	}
	*/
}

MshEditDialog::~MshEditDialog()
{
	delete _noDataDeleteBox;

	for (int i = 0; i < _labels.size(); i++)
	{
		delete _labels[i];
		delete _edits[i];
		delete _buttons[i];
	}
}

void MshEditDialog::accept()
{
/* TODO6
	int tabIndex = this->tabWidget->currentIndex();

	if (tabIndex >= 0)
	{
		MeshLib::CFEMesh* new_mesh = NULL;

		switch (tabIndex)
		{
		case 0:
		{
			const int nLayers = atoi(this->editNLayers->text().toStdString().c_str());
			const double thickness = strtod(replaceString(",", ".", this->editThickness->text().toStdString()).c_str(), 0);
			new_mesh = MshLayerMapper::CreateLayers(_msh, nLayers, thickness);
			break;
		}
		case 1:
		{
			new_mesh = new MeshLib::CFEMesh(*_msh);
			const size_t nLayers = _msh->getNumberOfMeshLayers();
			if (nLayers==0)
			{
				const std::string imgPath ( this->_edits[0]->text().toStdString() );
				if (!imgPath.empty())
					MshLayerMapper::LayerMapping(new_mesh, imgPath, nLayers, 0, _noDataDeleteBox->isChecked());
			}
			else
			{
				for (size_t i = 1; i <= nLayers+1; i++)
				{
					const std::string imgPath ( this->_edits[i]->text().toStdString() );
					if (!imgPath.empty())
					{
						int result = MshLayerMapper::LayerMapping(new_mesh, imgPath, nLayers, i-1, _noDataDeleteBox->isChecked());
						if (result==0) break;
					}
				}
			}
			if (nLayers>0 && this->_edits[0]->text().length()>0)
			{
				MeshLib::CFEMesh* final_mesh = MshLayerMapper::blendLayersWithSurface(new_mesh, nLayers, this->_edits[0]->text().toStdString());
				delete new_mesh;
				new_mesh = final_mesh;
			}
			break;
		}
		default:
			std::cout <<
			"Error in MshEditDialog::accept() - No instructions found for selected tab..."
			          << std::endl;
		}
		if (new_mesh)
		{
			std::string mshname("NewMesh");
			emit mshEditFinished(new_mesh, mshname);
		}
		else
			OGSError::box("Error creating mesh");
	}
	else
		std::cout << "Error in MshEditDialog::accept() - No tab selected... " << std::endl;
	this->done(QDialog::Accepted);
*/
}

void MshEditDialog::reject()
{
	this->done(QDialog::Rejected);
}

void MshEditDialog::getFileName()
{
	QPushButton* button = dynamic_cast<QPushButton*>(this->sender());
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString filename = QFileDialog::getOpenFileName(this,
	                                                "Select raster file to open",
	                                                settings.value(
	                                                        "lastOpenedFileDirectory").toString(
	                                                        ),
	                                                "ASCII raster files (*.asc);;All files (* *.*)");
	_fileButtonMap[button]->setText(filename);
	QDir dir = QDir(filename);
	settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
}
