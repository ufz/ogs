/**
 * \file DiagramPrefsDialog.h
 * KR Initial implementation
 */

#ifndef DIAGRAMPREFSDIALOG_H
#define DIAGRAMPREFSDIALOG_H

#include <QtGui/QMainWindow>
#include "ui_DiagramPrefs.h"
#include "DatabaseConnection.h"
#include "Station.h"

class DiagramList;

/**
 * \brief A dialog that allows for setting preferences for a requested diagram.
 *
 * A dialog that allows for setting preferences for a requested diagram. Note: In the current version
 * this dialog only works when requesting data from a database. Visualisation of data from an ASCII-file
 * is still possible using the "Load File"-button but setting the preferences will not work (i.e. it is
 * only possible to visualise all the data in the file with default preferences.
 */
class DiagramPrefsDialog : public QDialog, private Ui_DiagramPrefs
{
	Q_OBJECT

public:
	DiagramPrefsDialog(GEOLIB::Station* stn, QString listName, DatabaseConnection* db, QDialog* parent = 0);
	~DiagramPrefsDialog(void);


private:
	int loadFile(const QString &filename);
	int loadList(const std::vector< std::pair<QDateTime, float> > &coords);

	std::vector<DiagramList*> _list;
	DatabaseConnection* _db;
	int _listID;
	int _stationID;


private slots:
	void accept();
	void reject();
	void on_loadFileButton_clicked();

signals:


};

#endif //DIAGRAMPREFSDIALOG_H
