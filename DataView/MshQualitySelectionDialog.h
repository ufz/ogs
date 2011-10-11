/**
 * \file MshQualitySelectionDialog.h
 * 2011/03/16 KR Initial implementation
 */

#ifndef MSHQUALITYSELECTIONDIALOG_H
#define MSHQUALITYSELECTIONDIALOG_H

#include "MSHEnums.h"
#include "ui_MshQualitySelection.h"

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
