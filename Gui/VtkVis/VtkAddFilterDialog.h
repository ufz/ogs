/**
 * \file VtkAddFilterDialog.h
 * 23/2/2010 LB Initial implementation
 *
 */

#ifndef VTKADDFILTERDIALOG_H
#define VTKADDFILTERDIALOG_H

// ** INCLUDES **
#include "ui_VtkAddFilterDialogBase.h"

class VtkVisPipeline;
class QModelIndex;
class QRadioButton;

/**
 * \brief Dialog for selecting a filter to be applied to a VtkPipelineItem.
 * The dialog lets you select filters defined in VtkOGSFilter that have been registered as OGSFilterInfo - objects.
 */
class VtkAddFilterDialog : public QDialog, public Ui_VtkAddFilterDialogBase
{
	Q_OBJECT

public:
	/// Constructor
	VtkAddFilterDialog(VtkVisPipeline* pipeline, QModelIndex parentIndex, QDialog* parent = 0);

public slots:
	void on_buttonBox_accepted();

protected slots:
	void on_filterListWidget_currentRowChanged(int currentRow);

private:
	VtkVisPipeline* _pipeline;
	QModelIndex _parentIndex;
};

#endif // VTKADDFILTERDIALOG_H
