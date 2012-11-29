/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file VisPrefsDialog.h
 *
 * Created on 2010-06-14 by Karsten Rink
 */

#ifndef VISPREFSDIALOG_H
#define VISPREFSDIALOG_H

#include "Point.h"
#include "ui_VisPrefs.h"
#include <QDialog>

class VtkVisPipeline;
class VisualizationWidget;

/**
 * \brief A dialog window for settung up a visualisation preferences
 */
class VisPrefsDialog : public QDialog, private Ui_VisPrefsDialog
{
	Q_OBJECT

public:
	VisPrefsDialog(VtkVisPipeline* pipeline,
	               VisualizationWidget* widget,
	               QDialog* parent = NULL);

protected slots:
	/// @brief Sets the background colour.
	void on_bgColorButton_colorPicked(QColor color);

	/// @brief Adds a light above the scene.
	void on_lightAboveBox_clicked();

	/// @brief Adds a light below the scene.
	void on_lightBelowBox_clicked();

	/// @brief Sets the given superelevation on all vis pipeline source objects
	void on_superelevationPushButton_pressed();

	/// @brief
	void on_loadShowAllCheckBox_stateChanged(int state);

	/// @brief Culls backfacing rendering primitives on all actors.
	void on_cullBackfacesCheckBox_stateChanged(int state);

private:
	VtkVisPipeline* _vtkVisPipeline;
	VisualizationWidget* _visWidget;
	GeoLib::Point _above;
	GeoLib::Point _below;
};

#endif //VISPREFSDIALOG_H
