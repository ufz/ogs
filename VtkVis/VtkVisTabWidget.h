/**
 * \file VtkVisTabWidget.h
 * 18/2/2010 LB Initial implementation
 *
 */


#ifndef VTKVISTABWIDGET_H
#define VTKVISTABWIDGET_H

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

	void on_diffuseColorPickerButton_colorPicked(QColor color);
	void on_visibleEdgesCheckBox_stateChanged(int state);
	void on_edgeColorPickerButton_colorPicked(QColor color);
	void on_opacitySlider_sliderMoved(int value);
	void on_scaleZ_textChanged(const QString &text);

	void SetActiveAttributeOnItem(int idx);

private:
	void addColorTable();

	/// Reads the algorithm properties of the given pipeline item and builds a dialog for adjusting these properties in the GUI.
	void buildProportiesDialog(VtkVisPipelineItem* item);

	/// Reads the scalar arrays of the given vtk-object and constructs content for the scalar array selection box.
	void buildScalarArrayComboBox(vtkAlgorithm* algorithm);
	VtkVisPipelineItem* _item;

signals:
	/// Is emitted when a property was changed.
	void requestViewUpdate();

};

#endif // VTKVISTABWIDGET_H
