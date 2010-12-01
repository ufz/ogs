/**
 * \file DataViewWidget.cpp
 * 11/2/2010 LB Initial implementation
 * 
 * Implementation of DataViewWidget
 */

// ** INCLUDES **
#include "DataViewWidget.h"

#include "Model.h"

DataViewWidget::DataViewWidget( QWidget* parent /*= 0*/ )
: QWidget(parent)
{
	setupUi(this);

	dataView->setSelectionMode(QAbstractItemView::ExtendedSelection);
	dataView->setSelectionBehavior(QAbstractItemView::SelectRows);

	connect(this->modelSelectComboBox, SIGNAL(currentIndexChanged(int)),
		this, SLOT(changeActiveModel(int)));

	connect(this->removeAllPushButton, SIGNAL(clicked()), this, SLOT(emitRequestModelClear()));
}

void DataViewWidget::addModel( Model* model )
{
	this->modelSelectComboBox->addItem(model->name());
	_models.push_back(model);

	this->modelSelectComboBox->setCurrentIndex(-1); // hack for the first inserted model
	this->modelSelectComboBox->setCurrentIndex(_models.size() - 1);
}

void DataViewWidget::removeModel( Model* model )
{
	int index = 0;
	for (std::vector<Model*>::iterator it = _models.begin();
		it != _models.end(); ++it)
	{
		if (*it == model)
		{
			if (dataView->model() == model)
				dataView->setModel(NULL);

			_models.erase(it);
			modelSelectComboBox->removeItem(index);
			return;
		}
		index++;
	}
}

void DataViewWidget::changeActiveModel( int index )
{
	size_t numModels = _models.size() ;
	if (numModels == 0)
		return;
	if ((unsigned int)index > _models.size() - 1)
		return;

	dataView->setModel(_models[index]);
}

void DataViewWidget::emitRequestModelClear()
{
	std::string modelName = this->modelSelectComboBox->currentText().toStdString();
	emit requestModelClear(modelName);
}
