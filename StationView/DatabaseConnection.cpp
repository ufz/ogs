/**
 * \file DatabaseConnection.cpp
 * KR Initial implementation
 */

#include "DatabaseConnection.h"
#include <QSqlError>
#include <QSqlQuery>
#include <QSqlQueryModel>
#include <QVariant>
#include <QSettings>
#include "DateTools.h"
#include "OGSError.h"
#include "QueryResultsDialog.h"
#include "StringTools.h"

#include <iostream>
#include <fstream>


/// The OGS5-Connection to a database
DatabaseConnection::DatabaseConnection(GEOLIB::GEOObjects* geoObjects, QObject* parent) : QObject(parent), _geoObjects(geoObjects)
{
}

/// Deconstructor for the connection object
DatabaseConnection::~DatabaseConnection()
{
	_db.removeDatabase("QOCI");
}

/**
 * Initialising and testing the default database-connection.
 * \return 1 if the connection has been established, 0 (and an error message) otherwise.
 */
int DatabaseConnection::dbConnect()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString protocol = settings.value("DBProtocol", "").toString();
	QString hostname = settings.value("DBHost",     "").toString();
	QString dbname   = settings.value("DBName",     "").toString();
	QString user     = settings.value("DBUser",     "").toString();
	QString pass     = settings.value("DBPass",     "").toString();

	//default connection
	if (protocol.isEmpty() || hostname.isEmpty() || dbname.isEmpty() || user.isEmpty() || pass.isEmpty())
	{
		protocol = "QOCI";
		hostname = "cora1-vip.leipzig.ufz.de";
		dbname   = "ORACLE1";
		user     = "GEOSYS";
		pass     = "project09";
	}

	return dbConnect(protocol, hostname, dbname, user, pass);
}

/**
 * Initialising and testing a specific database-connection.
 * \param protocol The database connection driver
 * \param hostname The connection's host name
 * \param dbname The connection's database name
 * \param user The user name for connecting to the database
 * \param pass The password for the specified user name
 * \return 1 if the connection has been established, 0 (and an error message) otherwise.
 */
int DatabaseConnection::dbConnect(QString protocol, QString hostname, QString dbname, QString user, QString pass)
{
	QSqlDatabase::removeDatabase(_db.connectionName());

	_db = QSqlDatabase::addDatabase(protocol);
	_db.setHostName(hostname);
	_db.setDatabaseName(dbname);
	_db.setUserName(user);
	_db.setPassword(pass);

	return test(false);
}

/**
 * Setting up a new database-connection as default connection.
 * \param protocol The database connection driver
 * \param hostname The connection's host name
 * \param dbname The connection's database name
 * \param user The user name for connecting to the database
 * \param pass The password for the specified user name
 * \return 1 if the connection has been established, 0 (and an error message) otherwise.
 */
int DatabaseConnection::setConnection(QString protocol, QString hostname, QString dbname, QString user, QString pass)
{
	QSettings settings("UFZ", "OpenGeoSys-5");

	settings.setValue("DBProtocol", protocol);
	settings.setValue("DBHost",     hostname);
	settings.setValue("DBName",     dbname);
	settings.setValue("DBUser",     user);
	settings.setValue("DBPass",     pass);

	dbConnect(protocol, hostname, dbname, user, pass);
	disconnect(this, SLOT(setConnection(QString, QString, QString, QString, QString)));
	return test(true);
}

/**
 * Tests the current database connection.
 * \return 1 if the connection has been established, 0 (and an error message) otherwise.
 */
int DatabaseConnection::test(bool showBox)
{
	if (_db.open())
	{
		if (showBox) OGSError::box("Database connection established.");
		else std::cout << "Database connection established...\n";
		_db.close();
		return 1;
	}
	else
	{
		std::cout << "Database connection failed...\n";
		databaseError();
		OGSError::box("Could not connect to database.");
		return 0;
	}
}


/// Outputting an error related to the database connection
void DatabaseConnection::databaseError()
{
	QSqlError error = _db.lastError();
	if (error.isValid())
	{
		std::cout << (error.databaseText()).toStdString();
		std::cout << (error.driverText()).toStdString() << "\n\n";
	}
}

/// Opening a dialog containing all the available station lists
void DatabaseConnection::getListSelection()
{
	if (_db.open())
	{
		QSqlQueryModel *qModel = new QSqlQueryModel();
		qModel->setQuery("select l.listid, c.catname, l.listname from lists l, categories c where c.catid=l.catid order by c.catname");

		QueryResultsDialog* dbView = new QueryResultsDialog();
		connect(dbView, SIGNAL(listSelected(int)), this, SLOT(loadStationList(int)));

		dbView->setView(qModel);
		dbView->exec();

		_db.close();

	}
	else
		databaseError();
}

/// Tests if the current database connection is valid.
bool DatabaseConnection::isConnected()
{
	return _db.isValid();
}

/**
 * Loads a list of stations with a random colour
 * \param listID The ID of the list that is requested
 * \return 1 if there were no errors, 0 and an error message otherwise.
 */
int DatabaseConnection::loadStationList(const int &listID)
{

	return loadStationList(listID, GEOLIB::getRandomColor());
}

/**
 * Loads a list of stations with a specified colour
 * \param listID The ID of the list that is requested
 * \param color The colour which will be assigned to the objects in the list
 * \return 1 if there were no errors, 0 and an error message otherwise.
 */
int DatabaseConnection::loadStationList(int listID, const GEOLIB::Color* const color)
{
	if (_db.open())
	{
		QSqlQuery query, stnQuery;
		QString stationName;
		std::vector<GEOLIB::Point*> *stations = new std::vector<GEOLIB::Point*>;

		query.exec("select stationtype from geosysstationtypes where listid=" + QString::number(listID));

		if (query.next())
		{
			GEOLIB::Station::StationType type = static_cast<GEOLIB::Station::StationType>(query.value(0).toInt());
			if (type == GEOLIB::Station::BOREHOLE)
				query.exec("select c.catname, l.listname from lists l, categories c, boreholes b where c.catid=l.catid and l.listid=b.listid and l.listid=" + QString::number(listID));
			else
				query.exec("select c.catname, l.listname from lists l, categories c where c.catid=l.catid and l.listid=" + QString::number(listID));

			if (query.next())
			{
				QString listName = (query.value(0)).toString() + " (" + (query.value(1)).toString() + ")";

				if (type == GEOLIB::Station::BOREHOLE)
					stnQuery.exec("select s.stationid, s.name, s.x, s.y, s.z, b.bdepth, to_char(b.bdate, 'YYYY-MM-DD') from stations s, boreholes b where s.listid=b.listid and s.stationid=b.stationid and s.listid=" + QString::number(listID) + " order by stationid");
				else
					stnQuery.exec("select stationid, name, x, y, z from stations where listid=" + QString::number(listID) + " order by stationid");

				while (stnQuery.next())
				{
					stationName = stnQuery.value(1).toString();
					if (stationName.isEmpty()) stationName = "Station" + stnQuery.value(0).toString();

					GEOLIB::Station* newStation;
					if (type == GEOLIB::Station::BOREHOLE) newStation = GEOLIB::StationBorehole::createStation(stationName.toStdString(), stnQuery.value(2).toDouble(), stnQuery.value(3).toDouble(), stnQuery.value(4).toDouble(), stnQuery.value(5).toDouble(), stnQuery.value(6).toString().toStdString());
					else newStation = GEOLIB::Station::createStation(stationName.toStdString(), stnQuery.value(2).toDouble(), stnQuery.value(3).toDouble(), stnQuery.value(4).toDouble());
					newStation->setColor(color);
					stations->push_back(newStation);
				}

				if (type == GEOLIB::Station::BOREHOLE)
					//addStratigraphy(listID, _geoObjects->getStationVec(listName.toStdString()));
					addStratigraphy(listID, stations);

				std::string temp_name (listName.toStdString());
				if (!stations->empty()) _geoObjects->addStationVec(stations, temp_name , color);

				_db.close();
				//emit listLoaded(listName);

				return 1;
			}
		}
		else
		{
			std::cout << "DatabaseConnection::loadList() - No database entry found for the selected key." << std::endl;
			_db.close();
		}
	}
	else
		databaseError();

	return 0;
}

/**
 * Loads additional stratigraphy-data if the loaded station list consists of boreholes
 * \param listID The ID of the list that is requested
 * \param stations List of station objects for which stratigraphy data will be provided
 * \return 1 if there were no errors, 0 and an error message otherwise.
 */
int DatabaseConnection::addStratigraphy(int listID, std::vector<GEOLIB::Point*> *stations)
{
	if (_db.open())
	{
		QSqlQuery strat;

		size_t size = stations->size();
		for (size_t i=0; i<size; i++)
		{
			int count = 1;
			GEOLIB::StationBorehole* newStation = static_cast<GEOLIB::StationBorehole*>((*stations)[i]);
			strat.exec("select s.layerid, s.thickness, s.strat from stations t, stratigraphies s where s.listid=" + QString::number(listID) + " and s.listid=t.listid and s.stationid=t.stationid and t.name='" + QString::fromStdString(static_cast<GEOLIB::Station*>((*stations)[i])->getName()) + "' order by layerid");

			while (strat.next())
			{
				if (count==strat.value(0))
				{
					newStation->addSoilLayer (strat.value(1).toDouble(), (strat.value(2).toString()).toStdString());
					//newStation->type = Station::BOREHOLE;
				}
				else
				{
					std::cout << "DatabaseConnection::addStratigraphy - Station " << static_cast<GEOLIB::Station*>((*stations)[i])->getName() << ": Stratigraphy incomplete...\n";
				}
				count++;
			}
			(*stations)[i] = newStation;
		}
	}
	else
	{
		std::cout << "Database error" << std::endl;
		return 0;
	}

	return 1;
}

/**
 * Returns the list ID for the station at a given postition
 * \param list The name of the list
 * \param x The x-coordinate of a station within the list
 * \param y The y-coordinate of a station within the list
 * \return The list ID if there were no errors, -1 and an error message otherwise.
 */
int DatabaseConnection::getListID(const QString &list, const double &x, const double &y)
{
	if (_db.open())
	{
		QSqlQuery query;
		query.exec("select l.listid from lists l, categories c, stations s where l.catid = c.catid and l.listid = s.listid and s.x="
			       + QString::number(x,'f') + " and s.y=" + QString::number(y,'f') + " and c.catname='" + list.left(list.indexOf(" (")) + "'");

		if (query.next())
			return query.value(0).toInt();

		return -1;
	}

	return -1;
}

/**
 * Returns the station ID for the station at a given postition
 * \param listID The ID of the list the station belongs to
 * \param x The x-coordinate of the station
 * \param y The y-coordinate of the station
 * \return The station ID if there were no errors, -1 and an error message otherwise.
 */
int DatabaseConnection::getStationID(const int &listID, const double &x, const double &y)
{
	if (_db.open())
	{
		QSqlQuery query;
		query.exec("select stationid from stations where listid=" + QString::number(listID) + " and x=" + QString::number(x,'f') + " and y=" + QString::number(y,'f'));

		QString oldq = query.lastQuery();
		if (query.next())
			return query.value(0).toInt();

		return -1;
	}

	return -1;
}

/**
 * Returns all properties associated with a given list
 * Note: This function is currently not implemented correctly. Please leave it alone or contact me if you want to use it.
 * \param listID The ID of the queried list
 * \param propNames A vector in which the properties will be stored
 * \return 1 if there were no errors, 0 and an error message otherwise.
 */
int DatabaseConnection::getListProperties(const int &listID, std::vector<QString> &propNames)
{
	Q_UNUSED (listID);

	if (_db.open())
	{
		QSqlQuery query;
		query.exec("select column_name,data_type from all_tab_cols where table_name =\'STATIONS\'");
		while (query.next())
		{
			if ((query.value(0).toString()).compare("LISTID") != 0 && (query.value(0).toString()).compare("STATIONID") != 0 && (query.value(1).toString()).compare("VARCHAR2") != 0)
			{
				propNames.push_back(query.value(0).toString());
			}
		}
		return 1;
	}
	return 0;
}

/**
 * The minimum and maximum date for time series data associated with a given station
 * \param listID The ID of the list the station belongs to
 * \param stationID The ID of the station
 * \param startDate The value of the first date for that station found in the database
 * \param endDate The value of the last date for that station found in the database
 * \return 1 if there were no errors, 0 and an error message otherwise.
 */
int DatabaseConnection::getDateBounds(const int &listID, const int &stationID, QString &startDate, QString &endDate)
{
	startDate = "";
	endDate = "";
	if (_db.open())
	{
		QSqlQuery query;
		query.exec("select to_char(min(mdate), 'DD.MM.YYYY'), to_char(max(mdate), 'DD.MM.YYYY') from mvalues where listid=" + QString::number(listID) + " and stationID=" + QString::number(stationID));
		if (query.next())
		{
			startDate = query.value(0).toString();
			endDate = query.value(1).toString();
		}
		_db.close();
		return 1;
	}
	return 0;
}

/**
 * The minimum and maximum value for a given property associated with a given list
 * \param listID The ID of the list
 * \param prop The name of the property
 * \param min The smallest value of that property found in the database
 * \param max The largest value of that property found in the database
 * \return 1 if there were no errors, 0 and an error message otherwise.
 */
int DatabaseConnection::getPropertyBounds(const int &listID, const QString &prop, double &min, double &max)
{
	if (_db.open())
	{
		QSqlQuery query;
		query.exec("select min(" + prop + "),  max(" + prop + ") from stations where listid=" + QString::number(listID));
		if (query.next())
		{
			min = query.value(0).toDouble();
			max = query.value(1).toDouble();
		}
		_db.close();
		return 1;
	}
	return 0;
}

/**
 * Load time series data for a given station from the database
 * \param listID The ID of the list the station belongs to
 * \param stationID The ID of the station
 * \param startDate The start date for the requested time series data
 * \param endDate The end date for the requested time series data
 * \param xCoords x-coordinates of the retrieved data
 * \param yCoords y-coordinates of the retrieved data
 * \return 1 if there were no errors, 0 and an error message otherwise.
 */
int DatabaseConnection::loadValues(const int &listID, const int &stationID, const QDateTime &startDate, const QDateTime &endDate, std::vector< std::pair<QDateTime, float> > &values)
{
	if (startDate<endDate && _db.open())
	{
		QSqlQuery query;

		query.prepare("select to_char(mdate,'YYYY-MM-DD'),mvalue from mvalues where listid=:listID and stationid=:stationID "
					  "and mdate>=to_date(:startDate,'DD.MM.YYYY') and mdate<=to_date(:endDate,'DD.MM.YYYY') order by mdate");
		query.bindValue(":listID",		listID);
		query.bindValue(":stationID",	stationID);
		query.bindValue(":startDate",	startDate.toString("dd.MM.yyyy"));
		query.bindValue(":endDate",		endDate.toString("dd.MM.yyyy"));
		query.exec();

		while (query.next())
		{
			values.push_back( std::pair<QDateTime, float>(query.value(0).toDateTime(), static_cast<float>(query.value(1).toDouble())) );
		}

		_db.close();
		return 1;
	}

	_db.close();
	return 0;
}








/*********************************************
 * Inserting data into the database.         *
 * Be very careful what you do with any of   *
 * functions below.					         *
 * You might corrupt the database otherwise. *
 * --KR                                      *
 *********************************************/


/**
 * Inserts a new station list into the database
 * \param path		the	path to the file containing the data
 * \param listname	the name of the list given a certain category
 * \param catname	the category of stations (i.e. boreholes)
 * \param type		the OGS5 Stationtype
 */
int DatabaseConnection::addListToDB(std::string path, std::string listName, std::string catName, GEOLIB::Station::StationType type)
{
	QSqlQuery query;
	int listID, catID;
	std::string line;
	bool status=true, commit=true;

	if (_db.open())
	{
		std::ifstream in( path.c_str() );

		if (!in.is_open())
		{
			std::cout << "DatabaseConnection::addListToDB() - Could not open file..." << std::endl;
			return 0;
		}

		/* a name for the list in the first line of the file. this name is ignored here because
		 * the real list name is required as parameter to this method
		 */
		getline(in, line);
		if ((line.substr(0,1)).compare("!")==0)
			line.substr( 1, line.length()-1 );

		query.exec("select max(listid) from lists");
		if (query.next())
		{
			listID = query.value(0).toInt() + 1;

			query.exec("select catid from categories where catname='" + QString::fromStdString(catName) + "'");

			if (query.next())
				catID = query.value(0).toInt();
			else
			{
				query.exec("select max(catid) from categories");
				if (query.next())
					catID = query.value(0).toInt() + 1;
				query.prepare("insert into categories values (" + QString::number(catID) + ", '" + QString::fromStdString(catName) +"')");
				commitTransaction(query);
			}

			_db.transaction();
			query.exec("insert into lists values(" + QString::number(listID) + ", '" + QString::fromStdString(listName) + "', " + QString::number(catID) + ", 0)");
			if (type == GEOLIB::Station::BOREHOLE) query.exec("insert into geosysstationtypes values (" + QString::number(listID) + ", 2)");

			int stationID=1;

			/* read all stations */
			while ( getline(in, line) )
			{
				if (type == GEOLIB::Station::BOREHOLE) status = addBoreholeToDB(listID, stationID, line);
				else status = addStationToDB(listID, stationID, line);

				if (!status)
				{
					databaseError();
					commit=false;
				}
				stationID++;
			}

			if (commit) _db.commit();
			else _db.rollback();

			_db.close();
			in.close();

			return commit;
		}
		else
			std::cout << "Database error." << std::endl;

		_db.close();
	}
	else
		databaseError();

	return 0;
}

/**
 * Inserts a new station into the database (this is to be called from addListToDB())
 * \param listID	the ID of the list the station belongs to
 * \param stationID	the ID of the station
 * \param line		a line of text containing all the necessary information for the station (typically from a textfile)
 */
bool DatabaseConnection::addStationToDB(int listID, int stationID, std::string line)
{
	QSqlQuery query;
	GEOLIB::Station* station = GEOLIB::Station::createStation(line);
	query.prepare("insert into stations values(:listid, :stationid, :stationname, :x, :y, :z)");
	query.bindValue(":listid", listID);
	query.bindValue(":stationid", stationID);
	query.bindValue(":stationname", QString::fromStdString(station->getName()));
	query.bindValue(":x", (*station)[0]);
	query.bindValue(":y", (*station)[1]);
	query.bindValue(":z", (*station)[2]);
	return query.exec();
}

/**
 * Inserts a new borehole into the database (this is to be called from addListToDB())
 * This is kind of an expansion of addStationToDB() which uses the same parameters
 * \param listID	the ID of the list the station belongs to
 * \param stationID	the ID of the station
 * \param line		a line of text containing all the necessary information for the station (typically from a textfile)
 */
bool DatabaseConnection::addBoreholeToDB(int listID, int stationID, std::string line)
{
	QSqlQuery query;
	GEOLIB::StationBorehole* station = GEOLIB::StationBorehole::createStation(line);

	if (addStationToDB(listID, stationID, line))
	{
		query.prepare("insert into boreholes values (:listid, :stationid, :bdepth, to_date(:bdate, 'DD.MM.YYYY'))");
		query.bindValue(":listid", listID);
		query.bindValue(":stationid", stationID);
		query.bindValue(":bdepth", station->getDepth());
		QString sDate = QString::fromStdString(date2string(station->getDate()));
		query.bindValue(":bdate", sDate);
		return query.exec();
	}
	return false;
}

/**
 * Adds stratigraphic information for a given list of boreholes to the database
 * This method assumes that the given list of boreholes exists
 * \param path		the	path to the file containing the data
 * \param listID	the ID of the station list the stratigraphic data belongs to
 */
int DatabaseConnection::addStratigraphyToDB(std::string path, int listID)
{
	QSqlQuery query;
	int stationID;
	std::string line, stationName;

	if (_db.open())
	{
		std::ifstream in( path.c_str() );

		if (!in.is_open())
		{
			std::cout << "DatabaseConnection::addListToDB() - Could not open file..." << std::endl;
			return 0;
		}

		query.exec("select count(*) from lists where listid=" + QString::number(listID));
		if (query.next())
		{
			if (query.value(0).toInt() == 1)
			{

				_db.transaction();

				/* read all stations */
				while ( getline(in, line) )
				{
					std::list<std::string> fields = splitString(line, '\t');

					stationName = fields.front();
					fields.pop_front();
					query.exec("select stationid from stations where listid=" + QString::number(listID) + " and name='" + QString::fromStdString(stationName) + "'" );

					if (query.next() && fields.size() >= 3)
					{
						stationID = query.value(0).toInt();

						query.prepare("insert into stratigraphies values (:listid, :stationid, :layerid, :thickness, :strat, :petro)");
						query.bindValue(":listid", listID);
						query.bindValue(":stationid", stationID);
						query.bindValue(":layerid", atoi(fields.front().c_str()));
						fields.pop_front();
						query.bindValue(":thickness", strtod(replaceString(",", ".", fields.front()).c_str(),0));
						fields.pop_front();
						query.bindValue(":strat", QString::fromStdString(fields.front()));
						fields.pop_front();

						if (!fields.empty())
							query.bindValue(":petro", QString::fromStdString(fields.front()));
						else
							query.bindValue(":petro", "");

						query.exec();
					}
				}
				_db.commit();
			}
			in.close();
			_db.close();

			return 1;
		}
		else
			std::cout << "Database error." << std::endl;

		_db.close();
	}
	else
		databaseError();

	return 0;
}

/**
 * Adds time series information for a given station to the database
 * This method assumes that the given list and station-id exists
 * \param path		the	path to the file containing the data
 * \param listID	the ID of the list the stratigraphic data belongs to
 * \param stationID	the ID of the station the stratigraphic data belongs to
 */
int DatabaseConnection::addMeasuredValuesToDB(std::string path, int listID, int stationID)
{
	QSqlQuery query;
	std::string line;

	if (_db.open())
	{
		std::ifstream in( path.c_str() );

		if (!in.is_open())
		{
			std::cout << "DatabaseConnection::addMeasuredValuesToDB() - Could not open file..." << std::endl;
			return 0;
		}

		query.exec("select count(*) from stations where listid=" + QString::number(listID) + " and stationid=" + QString::number(stationID));
		if (query.next())
		{
			if (query.value(0).toInt() == 1)
			{

				_db.transaction();

				while ( getline(in, line) )
				{
					std::list<std::string> fields = splitString(line, '\t');

					QString a = query.lastQuery();

					query.prepare("insert into mvalues values (:listid, :stationid, to_date(:mdatetime,'DD.MM.YYYY'), :mvalue)");
					query.bindValue(":listid", listID);
					query.bindValue(":stationid", stationID);
					query.bindValue(":mdatetime", QString::fromStdString(fields.front()));
					fields.pop_front();
					query.bindValue(":mvalue", strtod(replaceString(",", ".", fields.front()).c_str(),0));
					fields.pop_front();
					query.exec();
				}
				_db.commit();
			}
			in.close();
			_db.close();

			return 1;
		}
		else
			std::cout << "Database error." << std::endl;

		_db.close();
	}
	else
		databaseError();

	return 0;
}



/**
 * Executing a query as an encapsulated transaction
 * \return True if everything is alright, false if there were errors.
 */
bool DatabaseConnection::commitTransaction(QSqlQuery query)
{
	_db.transaction();
	query.exec();
	bool r = _db.commit();
	return r;
}

