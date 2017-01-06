/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the DiagramPrefsDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DIAGRAMPREFSDIALOG_H
#define DIAGRAMPREFSDIALOG_H

#include "ui_DiagramPrefs.h"
#include <QMainWindow>


class DetailWindow;
class DiagramList;
class QCheckBox;

namespace GeoLib
{
class Station;
}

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
    /**
     * Opens a new dialog based on station and the list this station belongs to. If a database connection
     * is available, the program will try to find data associated with the station, otherwise data can be
     * loaded from a file.
     * \param stn The station object associated the diagram.
     * \param listName The station list the station belongs to.
     * \param parent The parent QDialog.
     */
    DiagramPrefsDialog(const GeoLib::Station* stn,
                       const QString &listName,
                       //DatabaseConnection* db,
                       QDialog* parent = 0);

    /**
     * Opens a new dialog and automatically reads data from the associated station object.
     * \param stn The station object associated the diagram.
     * \param parent The parent QDialog.
     */
    DiagramPrefsDialog(GeoLib::Station* stn, QDialog* parent = 0);

    /**
     * Opens a new dialog and automatically reads data from the specified file. The diagram is not associated
     * with any geometric object.
     * \param filename File containing data for the diagram(s) to be visualised.
     * \param[out] window Returns the created DetailWindow.
     * \param parent The parent QDialog.
     */
    DiagramPrefsDialog(const QString &filename,
                       DetailWindow* window = NULL,
                       QDialog* parent = 0);

    ~DiagramPrefsDialog(void);

private:
    /**
     * Creates checkboxes for every list of data values found. Per default all of these are checked, i.e. all
     * diagrams will be visualised. Any checkbox the user unchecks will result in the associated data not being
     * visualised.
     */
    void createVisibilityCheckboxes();

    /**
     * Loading data from a file
     * \param filename Name of the file containing the data
     * \return 1 if everything is okay, 0 and an error message if there were errors
     */
    int loadFile(const QString &filename);

    /**
     * Setting up the QDiagramList object were the time series data will be stored
     * \param coords List of coordinates.
     * \return 1 if everything is okay, 0 and an error message if there were errors
     */
    int loadList(const std::vector< std::pair<QDateTime, float> > &coords);

    std::vector<DiagramList*> _list;
    std::vector<QCheckBox*> _visability;
    int _listID;
    int _stationID;
    DetailWindow* _window;

private slots:
    /// Instructions if the OK-Button has been pressed.
    /// Note: Clicking the "Load from file"-button overrides the database input!
    void accept();

    /// Instructions if the Cancel-Button has been pressed.
    void reject();

    /// Instructions if the "Load File"-Button has been pressed.
    void on_loadFileButton_clicked();

signals:
};

#endif //DIAGRAMPREFSDIALOG_H
