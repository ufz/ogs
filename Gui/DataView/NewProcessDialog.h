/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file NewProcessDialog.h
 *
 * Created on 2011-11-17 by Karsten Rink
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
