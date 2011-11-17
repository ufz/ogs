/**
 * \file NewProcessDialog.h
 * 2011/11/17 KR Initial implementation
 */

#ifndef NEWPROCESSDIALOG_H
#define NEWPROCESSDIALOG_H

#include "ui_NewProcess.h"

#include <QDialog>

#include "ProjectData.h"

class ProcessInfo;

/**
 * \brief A dialog window for adding a new process in GUI
 */
class NewProcessDialog : public QDialog, private Ui_NewProcess
{
	Q_OBJECT

public:
	/// Constructor for creating a new FEM condition.
	NewProcessDialog(ProjectData &project, QDialog* parent = 0);

	~NewProcessDialog(void);

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
