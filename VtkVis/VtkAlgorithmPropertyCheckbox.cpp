/**
 * \file VtkAlgorithmPropertyCheckbox.cpp
 * 20/10/2010 LB Initial implementation
 *
 * Implementation of VtkAlgorithmPropertyCheckbox class
 */

// ** INCLUDES **
#include "VtkAlgorithmPropertyCheckbox.h"

#include "VtkAlgorithmProperties.h"

VtkAlgorithmPropertyCheckbox::VtkAlgorithmPropertyCheckbox(const bool value,
                                                           const QString& name,
                                                           VtkAlgorithmProperties* algProps,
                                                           QWidget* parent /*= 0*/ )
	: QCheckBox(parent), _name(name), _algProps(algProps)
{
	this->setChecked(value);
	connect(this, SIGNAL(stateChanged(int)), this, SLOT(setNewValue(int)));
}

VtkAlgorithmPropertyCheckbox::~VtkAlgorithmPropertyCheckbox()
{
}

void VtkAlgorithmPropertyCheckbox::setNewValue( int state )
{
	bool boolState = (bool)state;
	_algProps->SetUserProperty(_name, QVariant(boolState));
}
