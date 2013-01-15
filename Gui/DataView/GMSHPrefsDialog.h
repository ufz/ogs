/**
 * \file
 * \author Karsten Rink
 * \date   2010-01-21
 * \brief  Definition of the GMSHPrefsDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GMSHPREFSDIALOG_H
#define GMSHPREFSDIALOG_H

#include "ui_GMSHPrefs.h"
#include <QDialog>

class QStringListModel;

namespace GeoLib
{
class GEOObjects;
}

/**
 * \brief A dialog window for setting preferences for GMSH
 */
class GMSHPrefsDialog : public QDialog, private Ui_GMSHPrefs
{
	Q_OBJECT

public:
	GMSHPrefsDialog(const GeoLib::GEOObjects* geoObjects, QDialog* parent = 0);
	~GMSHPrefsDialog(void);

private:
	std::vector<std::string> getSelectedObjects(QStringList list);

	QStringListModel* _allGeo;
	QStringListModel* _selGeo;

private slots:
	void on_selectGeoButton_pressed();
	void on_deselectGeoButton_pressed();
	void on_radioAdaptive_toggled(bool isTrue);

	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

signals:
	void requestMeshing(std::vector<std::string> &, unsigned, double, double, double, bool);
};

#endif //GMSHPREFSDIALOG_H
