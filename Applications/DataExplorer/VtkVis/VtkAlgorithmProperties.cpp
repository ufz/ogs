/**
 * \file
 * \author Lars Bilke
 * \date   2012-02-14
 * \brief  Implementation of the VtkAlgorithmProperties class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"

#include <vtkProperty.h>
#include <vtkTexture.h>

#include "Applications/DataHolderLib/ColorLookupTable.h"
#include "Applications/FileIO/XmlIO/Qt/XmlLutReader.h"
#include "BaseLib/Logging.h"
#include "VtkColorLookupTable.h"

VtkAlgorithmProperties::VtkAlgorithmProperties(QObject* parent /*= nullptr*/)
    : QObject(parent)
{
    _property = vtkProperty::New();
    _texture = nullptr;
    _scalarVisibility = true;
    _algorithmUserProperties = new QMap<QString, QVariant>;
    _algorithmUserVectorProperties = new QMap<QString, QList<QVariant>>;
    _activeAttributeName = "";
    _removable = true;
}

VtkAlgorithmProperties::~VtkAlgorithmProperties()
{
    _property->Delete();
    if (_texture != nullptr)
    {
        _texture->Delete();
    }

    for (auto& row : _lut)
    {
        row.second->Delete();
    }
    delete _algorithmUserProperties;
    delete _algorithmUserVectorProperties;
}

vtkLookupTable* VtkAlgorithmProperties::GetLookupTable(
    const QString& array_name)
{
    auto it = _lut.find(array_name);
    if (it != _lut.end())
    {
        return it->second;
    }

    return nullptr;
}

void VtkAlgorithmProperties::RemoveLookupTable(const QString& array_name)
{
    auto it = _lut.find(array_name);
    if (it != _lut.end())
    {
        it->second->Delete();
        _lut.erase(it);
    }
}

void VtkAlgorithmProperties::SetLookUpTable(const QString& array_name,
                                            vtkLookupTable* lut)
{
    lut->Build();

    if (array_name.length() > 0)
    {
        this->RemoveLookupTable(array_name);
        _lut.insert(std::pair<QString, vtkLookupTable*>(array_name, lut));
        _activeAttributeName = array_name;
    }
}

void VtkAlgorithmProperties::SetLookUpTable(const QString& array_name,
                                            const QString& filename)
{
    DataHolderLib::ColorLookupTable lut;
    if (FileIO::XmlLutReader::readFromFile(filename, lut))
    {
        VtkColorLookupTable* colorLookupTable = VtkColorLookupTable::New();
        colorLookupTable->setLookupTable(lut);
        SetLookUpTable(array_name, colorLookupTable);
    }
    else
        ERR("Error reading color look-up table.");
}

void VtkAlgorithmProperties::SetScalarVisibility(bool on)
{
    _scalarVisibility = on;
    emit ScalarVisibilityChanged(on);
}

QVariant VtkAlgorithmProperties::GetUserProperty(QString name) const
{
    if (this->_algorithmUserProperties->contains(name))
    {
        return this->_algorithmUserProperties->value(name);
    }

    ERR("Not a valid property: {:s}", name.toStdString());
    return QVariant();
}

QList<QVariant> VtkAlgorithmProperties::GetUserVectorProperty(
    QString name) const
{
    if (this->_algorithmUserVectorProperties->contains(name))
    {
        return this->_algorithmUserVectorProperties->value(name);
    }

    ERR("Not a valid property: {:s}", name.toStdString());
    return QList<QVariant>();
}

void VtkAlgorithmProperties::SetActiveAttribute(QString name)
{
    if (name.contains("Solid Color") || name.contains("P-TextureCoordinates"))
    {
        SetScalarVisibility(false);
    }
    else
    {
        SetScalarVisibility(true);
    }
    _activeAttributeName = name;
}
