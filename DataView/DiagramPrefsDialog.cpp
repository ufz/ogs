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
 * \param listName
 * \param db The database connection were the diagram-related data can be found
 * \param parent The parent QDialog.
 */
DiagramPrefsDialog::DiagramPrefsDialog(GEOLIB::Station* stn, QString listName, DatabaseConnection* db, QDialog* parent) : QDialog(parent)
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
		DetailWindow* stationView = new DetailWindow(_list[0]);
		for (size_t i=1; i<_list.size(); i++)
			stationView->addList(_list[i]);
		stationView->setAttribute(Qt::WA_DeleteOnClose);
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
		QDateTime endDate = _list[0]->getStartDate().addSecs(_list[0]->maxXValue());
		toDateLine->setText(endDate.toString("dd.MM.yyyy"));//QString::number(_list[0]->maxXValue()));
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
