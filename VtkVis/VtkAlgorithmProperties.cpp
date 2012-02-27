/**
 * \file VtkAlgorithmProperties.cpp
 * 14/02/2012 LB Initial implementation
 * 
 * Implementation of the VtkAlgorithmProperties class
 */

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"

#include <vtkProperty.h>
#include <vtkTexture.h>

#include "VtkColorLookupTable.h"
#include "XmlIO/XmlLutReader.h"

VtkAlgorithmProperties::VtkAlgorithmProperties(QObject* parent /*= NULL*/)
	: QObject(parent)
{
	_property = vtkProperty::New();
	_texture  = NULL;
	_scalarVisibility = true;
	_algorithmUserProperties = new QMap<QString, QVariant>;
	_algorithmUserVectorProperties = new QMap<QString, QList<QVariant> >;
	_activeAttributeName = "";
}

VtkAlgorithmProperties::~VtkAlgorithmProperties()
{
	_property->Delete();
	if (_texture != NULL) _texture->Delete();

	for (std::map<QString, vtkLookupTable*>::iterator it = _lut.begin();
		it != _lut.end(); ++it)
		it->second->Delete();
	delete _algorithmUserProperties;
	delete _algorithmUserVectorProperties;
}

vtkLookupTable* VtkAlgorithmProperties::GetLookupTable(const QString& array_name)
{
	std::map<QString, vtkLookupTable*>::iterator it = _lut.find(array_name);
	if (it != _lut.end())
		return it->second;
	else
		return NULL;
}

void VtkAlgorithmProperties::SetLookUpTable(const QString &array_name, vtkLookupTable* lut)
{
	lut->Build();

	if (array_name.length() > 0)
	{
		std::map<QString, vtkLookupTable*>::iterator it = _lut.find(array_name);
		if (it != _lut.end())
		{
			it->second->Delete();
			_lut.erase(it);
		}
		_lut.insert( std::pair<QString, vtkLookupTable*>(array_name, lut) );
		_activeAttributeName = array_name;
	}
}

void VtkAlgorithmProperties::SetLookUpTable(const QString &array_name, const QString &filename)
{
	VtkColorLookupTable* colorLookupTable = XmlLutReader::readFromFile(filename);
	SetLookUpTable(array_name, colorLookupTable);
}

void VtkAlgorithmProperties::SetScalarVisibility(bool on)
{
	_scalarVisibility = on;
	emit ScalarVisibilityChanged(on);
}

QVariant VtkAlgorithmProperties::GetUserProperty(QString name) const
{
	if (this->_algorithmUserProperties->contains(name))
		return this->_algorithmUserProperties->value(name);
	else
	{
		std::cout << "Not a valid property: " << name.toStdString() << std::endl;
		return QVariant();
	}
}

QList<QVariant> VtkAlgorithmProperties::GetUserVectorProperty(QString name) const
{
	if (this->_algorithmUserVectorProperties->contains(name))
		return this->_algorithmUserVectorProperties->value(name);
	else
	{
		std::cout << "Not a valid property: " << name.toStdString() << std::endl;
		return QList<QVariant>();
	}
}
