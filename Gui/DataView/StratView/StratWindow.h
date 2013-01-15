/**
 * \file
 * \author Karsten Rink
 * \date   2010-03-16
 * \brief  Definition of the StratWindow class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef STRATWINDOW_H
#define STRATWINDOW_H

#include "ui_StratWindow.h"
#include <QWidget>

namespace GeoLib
{
class StationBorehole;
}

/**
 * \brief Creates a window to visualise the stratigraphy of a borehole.
 */
class StratWindow : public QWidget, public Ui_StratWindow
{
	Q_OBJECT

public:
	/**
	 * Constructor
	 * \param station The borehole object to be visualised.
	 * \param stratColors A color map.
	 * \param parent The parent QWidget.
	 */
	StratWindow(GeoLib::StationBorehole* station,
	            std::map<std::string, GeoLib::Color*>* stratColors = NULL,
	            QWidget* parent = 0);
	~StratWindow(void) { this->destroy(); }

private:
	/// Automatically resize window based on the measurements of the borehole.
	void resizeWindow();

private slots:
	void on_closeButton_clicked();
};

#endif //STRATWINDOW_H
