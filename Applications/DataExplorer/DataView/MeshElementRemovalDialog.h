/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-16
 * \brief  Definition of the MeshElementRemovalDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHELEMENTREMOVALDIALOG_H
#define MESHELEMENTREMOVALDIALOG_H

#include "ui_MeshElementRemoval.h"
#include <QDialog>

#ifndef Q_MOC_RUN  // See: https://bugreports.qt-project.org/browse/QTBUG-22829
#include "Applications/ApplicationsLib/ProjectData.h"
#endif

class Node;

/**
 * \brief A dialog window for settung up a database connection
 */
class MeshElementRemovalDialog : public QDialog, private Ui_MeshElementRemoval
{
	Q_OBJECT

public:
	MeshElementRemovalDialog(const ProjectData &project, QDialog* parent = 0);
	~MeshElementRemovalDialog(void);

private slots:
	void on_boundingBoxCheckBox_toggled(bool is_checked);
	void on_elementTypeCheckBox_toggled(bool is_checked);
	void on_materialIDCheckBox_toggled(bool is_checked);
	void on_meshNameComboBox_currentIndexChanged(int idx);
	void accept();
	void reject();

private:
	const ProjectData& _project;
	unsigned _currentIndex, _aabbIndex, _matIDIndex;

signals:
	void meshAdded(MeshLib::Mesh* mesh);
};

#endif //MESHELEMENTREMOVALDIALOG_H
