/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-23
 * \brief  Definition of the VtkAddFilterDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
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
    VtkAddFilterDialog(VtkVisPipeline &pipeline,
                       QModelIndex parentIndex,
                       QDialog* parent = nullptr);

public slots:
    void on_buttonBox_accepted();

protected slots:
    void on_filterListWidget_currentRowChanged(int currentRow);

private:
    VtkVisPipeline& _pipeline;
    QModelIndex _parentIndex;
};

#endif // VTKADDFILTERDIALOG_H
