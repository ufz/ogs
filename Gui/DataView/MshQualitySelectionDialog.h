/**
 * \file
 * \author Karsten Rink
 * \date   2011-03-16
 * \brief  Definition of the MshQualitySelectionDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MSHQUALITYSELECTIONDIALOG_H
#define MSHQUALITYSELECTIONDIALOG_H

#include "MshEnums.h"
#include "ui_MshQualitySelection.h"
#include <QDialog>

class VtkMeshSource;

/**
 * \brief A dialog window for settung up a database connection
 */
class MshQualitySelectionDialog : public QDialog, private Ui_MshQualitySelection
{
	Q_OBJECT

public:
	MshQualitySelectionDialog(VtkMeshSource* msh, QDialog* parent = 0);
	~MshQualitySelectionDialog(void);

private:
	VtkMeshSource* _msh;

private slots:
	void accept();
	void reject();

signals:
	void measureSelected(VtkMeshSource*, MshQualityType::type);
};

#endif //MSHQUALITYSELECTIONDIALOG_H
