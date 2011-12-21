/**
 * \file MshEditDialog.cpp
 * 2010/11/09 KR Initial implementation
 */

#include "MshEditDialog.h"
#include "OGSError.h"
#include "StringTools.h"
#include "msh_mesh.h"

#include <QCheckBox>
#include <QFileDialog>
#include <QPushButton>
#include <QSettings>

MshEditDialog::MshEditDialog(const MeshLib::CFEMesh* mesh, QDialog* parent)
	: QDialog(parent), _msh(mesh), _noDataDeleteBox(NULL)
{
	setupUi(this);

	this->gridLayoutLayerMapping->setMargin(5);
	this->gridLayoutLayerMapping->setColumnMinimumWidth(2,10);
	this->gridLayoutLayerMapping->setColumnStretch(0, 80);
	this->gridLayoutLayerMapping->setColumnStretch(1, 200);
	this->gridLayoutLayerMapping->setColumnStretch(2, 10);

	size_t nLayers = mesh->getNumberOfMeshLayers();
	if (nLayers == 0)
		nLayers = 1;   // adapt to old files where 2D meshes officially had "0" layers which makes no sense

	for (size_t i = 0; i < nLayers; i++)
	{
		QString text = (i) ? "Layer" + QString::number(i) : "Surface";
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
	}

	_noDataDeleteBox = new QCheckBox("Remove mesh nodes at NoData values");
	_noDataDeleteBox->setChecked(false);
	_noDataDeleteBox->setEnabled(false);
	if (nLayers == 1)
	{
		_noDataDeleteBox->setEnabled(true);
		this->gridLayoutLayerMapping->addWidget(_noDataDeleteBox, 2, 1);
	}
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
	int tabIndex = this->tabWidget->currentIndex();

	if (tabIndex >= 0)
	{
		MeshLib::CFEMesh* new_mesh = NULL;

		switch (tabIndex)
		{
		case 0:
		{
			int nLayers = atoi(this->editNLayers->text().toStdString().c_str());
			double thickness =
			        strtod(replaceString(",", ".", this->editThickness->text().toStdString()).
			               c_str(), 0);

			new_mesh = MshLayerMapper::CreateLayers(_msh, nLayers, thickness);
			break;
		}
		case 1:
		{
			size_t nLayers = _msh->getNumberOfMeshLayers();
			if (nLayers == 0)
				nLayers = 1;  // adapt to old files where 2D meshes officially had "0" layers which makes no sense

			for (size_t i = 0; i < nLayers; i++)
			{
				std::string imgPath ( this->_edits[i]->text().toStdString() );
				if (!imgPath.empty())
					new_mesh = MshLayerMapper::LayerMapping(
					        _msh,
					        imgPath,
					        nLayers,
					        i,
					        _noDataDeleteBox->
					        isChecked());
			}
			//if (nLayers>1) MshLayerMapper::CheckLayerMapping(new_mesh, nLayers, 1); //TODO !!!
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
