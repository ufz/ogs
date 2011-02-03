/**
 * \file DataViewWidget.h
 * 11/2/2010 LB Initial implementation
 * 
 */


#ifndef DATAVIEWWIDGET_H
#define DATAVIEWWIDGET_H

// ** INCLUDES **
#include "ui_DataViewWidgetBase.h"

class Model;

/**
 * \brief Widget containing views for geometrical objects.
 */
class DataViewWidget : public QWidget, public Ui_DataViewWidgetBase
{
	Q_OBJECT

public:
	DataViewWidget(QWidget* parent = 0);
	
public slots:
	void addModel(Model* model);
	void removeModel(Model* model);

private slots:
	void changeActiveModel(int index);
	void emitRequestModelClear();

signals:
	void requestModelClear(const std::string &name);

private:
	std::vector<Model*> _models;

};

#endif // DATAVIEWWIDGET_H

