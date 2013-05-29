/**
 * \file
 * \author Karsten Rink
 * \date   2013-05-29
 * \brief  Definition of the MergeGeometriesDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MERGEGEOMETRIESDIALOG_H
#define MERGEGEOMETRIESDIALOG_H

#include "ui_MergeGeometries.h"
#include <QDialog>

class QStringListModel;

namespace GeoLib
{
class GEOObjects;
}

/**
 * \brief A dialog window for setting preferences for GMSH
 */
class MergeGeometriesDialog : public QDialog, private Ui_MergeGeometries
{
	Q_OBJECT

public:
	MergeGeometriesDialog(GeoLib::GEOObjects* geoObjects, QDialog* parent = 0);
	~MergeGeometriesDialog(void);

private:
	std::vector<std::string> getSelectedGeometries(QStringList list);

	GeoLib::GEOObjects* _geo_objects;
	QStringListModel* _allGeo;
	QStringListModel* _selGeo;

private slots:
	void on_selectGeoButton_pressed();
	void on_deselectGeoButton_pressed();

	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

signals:
	void requestMeshing(std::vector<std::string> &, unsigned, double, double, double, bool);
};

#endif //MERGEGEOMETRIESDIALOG_H
