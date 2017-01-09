/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-18
 * \brief  Definition of the VtkVisTabWidget class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#pragma once

// ** INCLUDES **
#include "ui_VtkVisTabWidgetBase.h"

class vtkAlgorithm;

/**
 * \brief Contains a QTreeView of the VtkVisPipeline and a properties
 * panel for adjusting vtkAlgorithms rendering and filter settings.
 */
class VtkVisTabWidget : public QWidget, public Ui_VtkVisTabWidgetBase
{
    Q_OBJECT

public:
    /// Constructor
    VtkVisTabWidget(QWidget* parent = 0);

protected slots:
    /// Updates the property panels to show informations on the given VtkVisPipelineItem.
    void setActiveItem(VtkVisPipelineItem* item);

    void on_arrayResetPushButton_clicked();
    void on_diffuseColorPickerButton_colorPicked(QColor color);
    void on_visibleEdgesCheckBox_stateChanged(int state);
    void on_edgeColorPickerButton_colorPicked(QColor color);
    void on_opacitySlider_sliderMoved(int value);
    void on_scaleZ_textChanged(const QString &text);
    void on_transX_textChanged(const QString &text)
    {
        Q_UNUSED(text);
        this->translateItem();
    }
    void on_transY_textChanged(const QString &text)
    {
        Q_UNUSED(text);
        this->translateItem();
    }
    void on_transZ_textChanged(const QString &text)
    {
        Q_UNUSED(text);
        this->translateItem();
    }

    void SetActiveAttributeOnItem(const QString &name);

private:
    /// Reads the algorithm properties of the given pipeline item and builds a dialog for adjusting these properties in the GUI.
    void buildProportiesDialog(VtkVisPipelineItem* item);

    /// Reads the scalar arrays of the given vtk-object and constructs content for the scalar array selection box.
    void buildScalarArrayComboBox(VtkVisPipelineItem* item);

    void translateItem();

    VtkVisPipelineItem* _item;

signals:
    /// Is emitted when a property was changed.
    void requestViewUpdate();
};
