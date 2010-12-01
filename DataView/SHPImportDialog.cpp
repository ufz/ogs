/**
 * \file SHPImportDialog.cpp
 * 25/01/2010 KR Initial implementation
 */

#include "OGSError.h"
#include "SHPImportDialog.h"
#include "SHPInterface.h"
#include "GEOModels.h"

#include <QLabel>
#include <QVBoxLayout>
#include <QFileInfo>
#include <QDialogButtonBox>
#include <QLineEdit>
#include <QRadioButton>


SHPImportDialog::SHPImportDialog(std::string filename, GEOModels* geoModels, QDialog* parent) 
: QDialog(parent), _buttonBox(NULL), _layout(NULL), _shpContentLabel(NULL), _nameLabel(NULL),
_listName(new QLineEdit()), _choice1(NULL), _choice2(NULL), _filename(filename), _fileType(0), 
_shpInterface(new SHPInterface(geoModels))
{
	setupDialog();
	show();
}


SHPImportDialog::~SHPImportDialog()
{
	delete _shpInterface;
	delete _buttonBox;
	delete _listName;
	delete _choice1;
	delete _choice2;
	delete _shpContentLabel;
	delete _nameLabel;
	delete _layout;
}


void SHPImportDialog::setupDialog()
{
	_layout = new QGridLayout(this);
	int shpType = 0, numberOfEntities = 0;
	QString type = "";

	setWindowTitle("Import SHP File");

	if (_shpInterface->readSHPInfo(_filename, shpType, numberOfEntities))
	{
		if ((shpType-1)%10 == 0) type="points";
		if ((shpType-3)%10 == 0) type="polylines";
		if ((shpType-5)%10 == 0) type="polygons";
		if ((shpType-8)%10 == 0) type="multipoints";
		if (shpType == 31) type="TIN elements";

		_shpContentLabel = new QLabel("The selected file contains " + QString::number(numberOfEntities) + " " + type, this);
		_nameLabel = new QLabel("List Name: ", this);

		QFileInfo fi(QString::fromStdString(_filename));
		_listName->setText(fi.baseName());

		if ((shpType-1)%10 == 0 && shpType!=31)	// Points
		{
				_choice1 = new QRadioButton("Read as GLI-Points");
				_choice2 = new QRadioButton("Read as Station Points");
				_choice1->toggle(); // default choice
				_layout->addWidget( _shpContentLabel );
				_layout->addWidget( _choice1 );
				_layout->addWidget( _choice2 );
				_layout->addWidget( _nameLabel );
				_layout->addWidget( _listName );
				_fileType = 1;
		}
		else if ((shpType-3)%10 == 0 || (shpType-5)%10 == 0)	// Polylines
		{
				_choice1 = new QRadioButton("Read as Polylines");
				_choice2 = new QRadioButton("Read as Polygons");
				if ((shpType-3)%10 == 0) _choice2->setDisabled(true); // disable polygon-choice if file contains only polylines
				_choice1->toggle(); // default choice
				_layout->addWidget( _shpContentLabel );
				_layout->addWidget( _choice1 );
				_layout->addWidget( _choice2 );
				_layout->addWidget( _nameLabel );
				_layout->addWidget( _listName );
				_fileType = 2;
		}
		else
		{
				_nameLabel->setText("This element type is currently not supported.");
				_layout->addWidget( _shpContentLabel );
				_layout->addWidget( _nameLabel );
		}

		_buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
		connect(_buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
		connect(_buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
		_layout->addWidget( _buttonBox );

		setLayout(_layout);
	}
	else
	{
		OGSError::box("Error reading shapefile!");
	}
}

void SHPImportDialog::accept()
{
	if (_listName->text().compare("")==0)
	{
		OGSError::box("Please insert a name for the data in this file.");
		return;
	}
	else
	{
		if (_fileType==1 && _choice1->isChecked()) _shpInterface->readSHPFile(_filename, SHPInterface::POINT,    _listName->text().toStdString());
		if (_fileType==1 && _choice2->isChecked()) _shpInterface->readSHPFile(_filename, SHPInterface::STATION,  _listName->text().toStdString());
		if (_fileType==2 && _choice1->isChecked()) _shpInterface->readSHPFile(_filename, SHPInterface::POLYLINE, _listName->text().toStdString());
		if (_fileType==2 && _choice2->isChecked()) _shpInterface->readSHPFile(_filename, SHPInterface::POLYGON,  _listName->text().toStdString());
		emit shpLoaded(_listName->text());
	}
	this->done(QDialog::Accepted);
}

void SHPImportDialog::reject()
{
	this->done(QDialog::Rejected);
}
