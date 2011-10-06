/**
 * \file VtkAlgorithmPropertyVectorEdit.cpp
 * 22/10/2010 LB Initial implementation
 *
 * Implementation of VtkAlgorithmPropertyVectorEdit class
 */

// ** INCLUDES **
#include "VtkAlgorithmPropertyVectorEdit.h"

#include "VtkAlgorithmProperties.h"

#include <QDoubleValidator>
#include <QHBoxLayout>
#include <QIntValidator>
#include <QLineEdit>
#include <QSize>

VtkAlgorithmPropertyVectorEdit::VtkAlgorithmPropertyVectorEdit( const QList<QString> contents,
                                                                const QString& name,
                                                                QVariant::Type type,
                                                                VtkAlgorithmProperties* algProps,
                                                                QWidget* parent /*= 0*/ )
	: QWidget(parent), _name(name), _algProps(algProps), _type(type)
{
	QHBoxLayout* layout = new QHBoxLayout;
	layout->setSpacing(3);
	layout->setContentsMargins(0, 0, 0, 0);

	foreach(QString content, contents)
	{
		QLineEdit* lineEdit = new QLineEdit(content, this);
		layout->addWidget(lineEdit);

		switch(_type)
		{
		case QVariant::Double:
			lineEdit->setValidator(new QDoubleValidator(this));
			break;

		case QVariant::Int:
			lineEdit->setValidator(new QIntValidator(this));

		default:
			break;
		}

		connect(lineEdit, SIGNAL(editingFinished()), this, SLOT(setNewValue()));
	}

	this->setLayout(layout);
	this->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
}

VtkAlgorithmPropertyVectorEdit::~VtkAlgorithmPropertyVectorEdit()
{
}

void VtkAlgorithmPropertyVectorEdit::setNewValue()
{
	QLayout* layout = this->layout();
	QList<QVariant> list;
	for (int i = 0; i < layout->count(); ++i)
	{
		QLineEdit* lineEdit = static_cast<QLineEdit*>(layout->itemAt(i)->widget());
		list.push_back(QVariant(lineEdit->text()));
	}

	_algProps->SetUserVectorProperty(_name, list);

	emit editingFinished();
}
