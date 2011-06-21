/**
 * \file DiagramPrefsDialog.cpp
 * KR Initial implementation
 */

#include "DiagramPrefsDialog.h"
#include "DatabaseConnection.h"
#include "DetailWindow.h"
#include "DiagramList.h"
#include "OGSError.h"
#include "Station.h"

#include <QCheckBox>
#include <QFileDialog>
#include <QMessageBox>


/**
 * Opens a new dialog.
 * \param stn The station object associated the diagram
 * \param listName
 * \param db The database connection were the diagram-related data can be found
 * \param parent The parent QDialog.
 */
DiagramPrefsDialog::DiagramPrefsDialog(GEOLIB::Station* stn, QString listName, DatabaseConnection* db, QDialog* parent) 
: QDialog(parent)
{
	setAttribute(Qt::WA_DeleteOnClose);

	_listID = -1; _stationID = -1;
	setupUi(this);
	stationNameLabel->setText(QString::fromStdString(stn->getName()));
	stationTypeLabel->setText(listName);

	if (db)
	{
		_db = db;
		_listID    = _db->getListID(listName, (*stn)[0], (*stn)[1]);
		_stationID = _db->getStationID(_listID, (*stn)[0], (*stn)[1]);
		if (_listID > 0 && _stationID > 0)
		{
			QString startDate, endDate;
			if (_db->getDateBounds(_listID, _stationID, startDate, endDate))
			{
				fromDateLine->setText(startDate);
				toDateLine->setText(endDate);
			}
		}
	}
}

DiagramPrefsDialog::DiagramPrefsDialog(const QString &filename, QDialog* parent)
: QDialog(parent)
{
	QFileInfo fi(filename);
	setupUi(this);
	stationNameLabel->setText(fi.baseName());
	stationTypeLabel->setText("");
	this->loadFile(filename);
}

DiagramPrefsDialog::~DiagramPrefsDialog()
{
	for (size_t i=0; i<_list.size(); i++)
		delete _list[i];
	this->destroy();
}

/// Instructions if the OK-Button has been pressed.
/// Note: Clicking the "Load from file"-button overrides the database input!
void DiagramPrefsDialog::accept()
{
	if ((fromDateLine->text().length()>0) && (toDateLine->text().length()>0) && (!_list.empty()))
	{
		if (_list[0]->size() == 0)	// data will be read from the database (if data has been loaded from file, size is already >0)
		{
			if (_listID > 0 && _stationID > 0)
			{
				std::vector< std::pair<QDateTime, float> > values;
				_db->loadValues(_listID, _stationID, QDateTime::fromString(fromDateLine->text(), "dd.MM.yyyy"), QDateTime::fromString(toDateLine->text(), "dd.MM.yyyy"), values);
				if (!loadList(values))
					OGSError::box("No data found.");
			}
		}
		
		// data has been loaded
		if (_list[0]->size()>0)	
		{
			DetailWindow* stationView = new DetailWindow();
			stationView->setAttribute(Qt::WA_DeleteOnClose);
			bool window_is_empty(true);
			for (size_t i=0; i<_list.size(); i++)
			{
				if (this->_visability[i]->isChecked())
				{
					stationView->addList(_list[i]);
					window_is_empty = false;
				}
			}

			if (!window_is_empty)
			{
				stationView->show();
				this->done(QDialog::Accepted);
			}
			else
			{
				delete stationView;
				OGSError::box("No dataset selected.");
			}
		}
		else 
		{
			OGSError::box("Invalid station data.");
			this->done(QDialog::Rejected);
		}
	}
	else
		OGSError::box("No data found...");
}

/// Instructions if the Cancel-Button has been pressed.
void DiagramPrefsDialog::reject()
{
	this->done(QDialog::Rejected);
}

/// Instructions if the "Load File"-Button has been pressed.
void DiagramPrefsDialog::on_loadFileButton_clicked()
{
	QString fileName = QFileDialog::getOpenFileName(this, "Select time series file to open", "","Time series files (*.stn *.txt)");
	if (!fileName.isEmpty())
		loadFile(fileName);
}

/**
 * Loading data from a file
 * \param filename Name of the file containing the data
 * return 1 if everything is okay, 0 and an error message if there were errors
 */
int DiagramPrefsDialog::loadFile(const QString &filename)
{
	if (DiagramList::readList(filename, _list))
	{
		for (size_t i=0; i<_list.size(); i++)
		{
			//_list[i]->setName(stationTypeLabel->text() + ": " + stationNameLabel->text());
			_list[i]->setXLabel("Time");
			//_list[i]->setYLabel("Water Level");
			_list[i]->setXUnit("day");
			//_list[i]->setYUnit("metres");
			_list[i]->setColor(QColor(Qt::red));
		}
		fromDateLine->setText(_list[0]->getStartDate().toString("dd.MM.yyyy"));//QString::number(_list[0]->minXValue()));
		QDateTime endDate = _list[0]->getStartDate().addSecs(static_cast<int>(_list[0]->maxXValue()));
		toDateLine->setText(endDate.toString("dd.MM.yyyy"));//QString::number(_list[0]->maxXValue()));
		this->createVisibilityCheckboxes();
		return 1;
	}

	OGSError::box("Error reading file.");
	return 0;
}

/**
 * Setting up the QDiagramList object were the time series data will be stored
 * \param coords List of coordinates.
 * return 1 if everything is okay, 0 and an error message if there were errors
 */
int DiagramPrefsDialog::loadList(const std::vector< std::pair<QDateTime, float> > &coords)
{
	if (!coords.empty())
	{
		DiagramList* l = new DiagramList;
		l->setName(stationTypeLabel->text() + ": " + stationNameLabel->text());
		l->setXLabel("Time");
		//l->setYLabel("Water Level");
		l->setXUnit("day");
		//l->setYUnit("metres");
		l->setColor(QColor(Qt::red));
		l->setList(coords);
		_list.push_back(l);
		return 1;
	}
	return 0;
}

void DiagramPrefsDialog::createVisibilityCheckboxes()
{
	for (size_t i=0; i<_list.size(); i++)
	{
		QCheckBox* box = new QCheckBox(_list[i]->getName());
		box->setChecked(true);
		this->CheckBoxLayout->addWidget(box);
		_visability.push_back(box);
	}
}