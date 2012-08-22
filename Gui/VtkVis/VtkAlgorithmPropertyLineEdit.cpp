/**
 * \file VtkAlgorithmPropertyLineEdit.cpp
 * 18/10/2010 LB Initial implementation
 *
 * Implementation of VtkAlgorithmPropertyLineEdit class
 */

// ** INCLUDES **
#include "VtkAlgorithmPropertyLineEdit.h"

#include "VtkAlgorithmProperties.h"

#include <QDoubleValidator>
#include <QIntValidator>

VtkAlgorithmPropertyLineEdit::VtkAlgorithmPropertyLineEdit(const QString& contents,
                                                           const QString& name,
                                                           QVariant::Type type,
                                                           VtkAlgorithmProperties* algProps,
                                                           QWidget* parent /*= 0*/)
	: QLineEdit(contents, parent), _name(name), _algProps(algProps), _type(type)
{
	switch(_type)
	{
	case QVariant::Double:
		this->setValidator(new QDoubleValidator(this));
		break;

	case QVariant::Int:
		this->setValidator(new QIntValidator(this));

	default:
		break;
	}

	connect(this, SIGNAL(editingFinished()), this, SLOT(setNewValue()));
}

VtkAlgorithmPropertyLineEdit::~VtkAlgorithmPropertyLineEdit()
{
}

void VtkAlgorithmPropertyLineEdit::setNewValue()
{
	QVariant value(this->text());
	if (value.convert(_type))
		_algProps->SetUserProperty(_name, value);
}
