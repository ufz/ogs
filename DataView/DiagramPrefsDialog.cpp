/**
 * \file DiagramPrefsDialog.cpp
 * KR Initial implementation
 */

#include <QFileDialog>
#include <QMessageBox>
#include "DiagramPrefsDialog.h"
#include "DetailWindow.h"
#include "DiagramList.h"
#include "OGSError.h"


/**
 * Opens a new dialog.
 * \param stn The station object associated the diagram
 * \param db The database connection were the diagram-related data can be found
 */
DiagramPrefsDialog::DiagramPrefsDialog(GEOLIB::Station* stn, QString listName, DatabaseConnection* db, QDialog* parent) : QDialog(parent)
{
	setAttribute(Qt::WA_DeleteOnClose);

	_listID = -1; _stationID = -1;
	_list = new DiagramList();
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

DiagramPrefsDialog::~DiagramPrefsDialog()
{
	delete _list;
	this->destroy();
}

/// Instructions if the OK-Button has been pressed.
/// Note: Clicking the "Load from file"-button overrides the database input!
void DiagramPrefsDialog::accept()
{
	if (_list->size() == 0)	// data will be read from the database (if data has been loaded from file, size is already >0)
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
	if (_list->size()>0)	
	{
		DetailWindow* stationView = new DetailWindow(_list);
		stationView->show();
	}
	else 
		OGSError::box("Invalid station.");
	this->done(QDialog::Accepted);
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
	_list->setName(stationTypeLabel->text() + ": " + stationNameLabel->text());
	_list->setXLabel("Time");
	//_list->setYLabel("Water Level");
	_list->setXUnit("day");
	//_list->setYUnit("metres");
	_list->setColor(QColor(Qt::red));

	if (_list->readList(filename))
	{
		fromDateLine->setText(QString::number(_list->minXValue()));
		toDateLine->setText(QString::number(_list->maxXValue()));
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
		_list->setName(stationTypeLabel->text() + ": " + stationNameLabel->text());
		_list->setXLabel("Time");
		//_list->setYLabel("Water Level");
		_list->setXUnit("day");
		//_list->setYUnit("metres");
		_list->setColor(QColor(Qt::red));
		_list->setList(coords);
		return 1;
	}
	return 0;
}
