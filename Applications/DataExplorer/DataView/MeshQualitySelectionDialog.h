/**
 * \file
 * \author Karsten Rink
 * \date   2011-03-16
 * \brief  Definition of the MshQualitySelectionDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MSHQUALITYSELECTIONDIALOG_H
#define MSHQUALITYSELECTIONDIALOG_H

#include "MeshEnums.h"
#include "ui_MeshQualitySelection.h"
#include <QDialog>

class VtkMeshSource;

/**
 * \brief A dialog for selecting a mesh quality metric
 */
class MshQualitySelectionDialog : public QDialog, private Ui_MeshQualitySelection
{
	Q_OBJECT

public:
	MshQualitySelectionDialog(QDialog* parent = 0);
	~MshQualitySelectionDialog(void);

	MeshQualityType getSelectedMetric() { return _metric; }

private:
	MeshQualityType _metric;

private slots:
	void accept();
	void reject();
};

#endif //MSHQUALITYSELECTIONDIALOG_H
