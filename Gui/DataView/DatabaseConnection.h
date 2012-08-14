/**
 * \file DatabaseConnection.h
 * KR Initial implementation
 */

#ifndef DATABASECONNECTION_H
#define DATABASECONNECTION_H

#include "GEOObjects.h"
#include "Station.h"
#include <QDateTime>
#include <QObject>
#include <QSqlDatabase>
#include <string>

/**
 * \brief Management of a database connection including the actual connection process, error handling and relevant queries.
 */
class DatabaseConnection : public QObject
{
	Q_OBJECT

public:
	DatabaseConnection(GEOLIB::GEOObjects* geoObjects, QObject* parent = 0);
	~DatabaseConnection();

	int dbConnect();

	int getDateBounds(const int &listID,
	                  const int &stationID,
	                  QString &startDate,
	                  QString &endDate);
	int getListID(const QString &list, const double &x, const double &y);
	int getListProperties(const int &listID, std::vector<QString> &propNames);
	int getPropertyBounds(const int &listID, const QString &prop, double &min, double &max);
	int getStationID(const int &listID, const double &x, const double &y);
	void getListSelection();
	bool isConnected();
	int loadValues(const int &listID,
	               const int &stationID,
	               const QDateTime &startDate,
	               const QDateTime &endDate,
	               std::vector< std::pair<QDateTime, float> > &values);
	int test(bool showBox);

public slots:
	int dbConnect(QString protocol,
	              QString hostname,
	              QString dbname,
	              QString user,
	              QString pass);
	int loadStationList(int listID);
	int setConnection(QString protocol,
	                  QString hostname,
	                  QString dbname,
	                  QString user,
	                  QString pass);

private:
	void databaseError();
	int addStratigraphy(int listID, std::vector<GEOLIB::Point*>* stations);
	bool commitTransaction(QSqlQuery query);

	//data-insert functions -- be careful!
	int addListToDB(std::string path,
	                std::string listName,
	                std::string catName,
	                GEOLIB::Station::StationType type);
	bool addStationToDB(int listID, int stationID, std::string line);
	bool addBoreholeToDB(int listID, int stationID, std::string line);
	int addStratigraphyToDB(std::string path, int listID);
	int addMeasuredValuesToDB(std::string path, int listID, int stationID);

	QSqlDatabase _db;
	GEOLIB::GEOObjects* _geoObjects;

signals:
	void listLoaded(QString listName);
};
#endif //DATABASECONNECTION_H
