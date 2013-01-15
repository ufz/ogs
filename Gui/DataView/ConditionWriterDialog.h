/**
 * \file
 * \author Karsten Rink
 * \date   2012-01-11
 * \brief  Definition of the ConditionWriterDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CONDITIONWRITERDIALOG_H
#define CONDITIONWRITERDIALOG_H

#include "ui_ConditionWriter.h"
#include <QDialog>

#include "GEOObjects.h"
#include "FEMCondition.h"

/**
 * \brief A dialog window for creating DIRECT boundary conditions from raster files
 */
class ConditionWriterDialog : public QDialog, private Ui_ConditionWriter
{
	Q_OBJECT

public:
	ConditionWriterDialog(const GeoLib::GEOObjects* geoObjects, QDialog* parent = 0);
	~ConditionWriterDialog(void);

private slots:
	void on_fileNameButton_pressed();

	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

signals:
	void saveFEMConditionsRequested(const QString&, const FEMCondition::CondType, const QString&);

};

#endif //CONDITIONWRITERDIALOG_H
