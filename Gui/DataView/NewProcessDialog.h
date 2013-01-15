/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-17
 * \brief  Definition of the NewProcessDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NEWPROCESSDIALOG_H
#define NEWPROCESSDIALOG_H

#include <QDialog>

#include "ui_NewProcess.h"

class ProcessInfo;

/**
 * \brief A dialog window for adding a new process in GUI
 */
class NewProcessDialog : public QDialog, private Ui_NewProcess
{
	Q_OBJECT

public:
	/// Constructor for creating a new FEM condition.
	NewProcessDialog(QDialog* parent = 0);

	~NewProcessDialog(void) {};

private:
	void setupDialog();

private slots:
	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

signals:
	void addProcess(ProcessInfo*);

};

#endif //NEWPROCESSDIALOG_H
