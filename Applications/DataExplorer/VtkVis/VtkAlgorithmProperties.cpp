/**
 * \file
 * \author Lars Bilke
 * \date   2012-02-14
 * \brief  Implementation of the VtkAlgorithmProperties class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"

#include "BaseLib/Logging.h"

#include <vtkProperty.h>
#include <vtkTexture.h>

#include "Applications/DataHolderLib/ColorLookupTable.h"
#include "VtkColorLookupTable.h"
#include "XmlIO/Qt/XmlLutReader.h"

VtkAlgorithmProperties::VtkAlgorithmProperties(QObject* parent /*= nullptr*/)
    : QObject(parent)
{
    property_ = vtkProperty::New();
    texture_  = nullptr;
    scalarVisibility_ = true;
    algorithmUserProperties_ = new QMap<QString, QVariant>;
    algorithmUserVectorProperties_ = new QMap<QString, QList<QVariant> >;
    activeAttributeName_ = "";
    removable_ = true;
}

VtkAlgorithmProperties::~VtkAlgorithmProperties()
{
    property_->Delete();
    if (texture_ != nullptr)
    {
        texture_->Delete();
    }

    for (auto& row : lut_)
    {
        row.second->Delete();
    }
    delete algorithmUserProperties_;
    delete algorithmUserVectorProperties_;
}

vtkLookupTable* VtkAlgorithmProperties::GetLookupTable(const QString& array_name)
{
    auto it = lut_.find(array_name);
    if (it != lut_.end())
    {
        return it->second;
    }

    return nullptr;
}

void VtkAlgorithmProperties::RemoveLookupTable(const QString& array_name)
{
    auto it = lut_.find(array_name);
    if (it != lut_.end())
    {
        it->second->Delete();
        lut_.erase(it);
    }
}

void VtkAlgorithmProperties::SetLookUpTable(const QString &array_name, vtkLookupTable* lut)
{
    lut->Build();

    if (array_name.length() > 0)
    {
        this->RemoveLookupTable(array_name);
        lut_.insert( std::pair<QString, vtkLookupTable*>(array_name, lut) );
        activeAttributeName_ = array_name;
    }
}

void VtkAlgorithmProperties::SetLookUpTable(const QString &array_name, const QString &filename)
{
    DataHolderLib::ColorLookupTable lut;
    if (FileIO::XmlLutReader::readFromFile(filename, lut))
    {
        VtkColorLookupTable* colorLookupTable = VtkColorLookupTable::New();
        colorLookupTable->setLookupTable(lut);
        SetLookUpTable(array_name, colorLookupTable);
    }
    else
        ERR ("Error reading color look-up table.");
}

void VtkAlgorithmProperties::SetScalarVisibility(bool on)
{
    scalarVisibility_ = on;
    emit ScalarVisibilityChanged(on);
}

QVariant VtkAlgorithmProperties::GetUserProperty(QString name) const
{
    if (this->algorithmUserProperties_->contains(name))
    {
        return this->algorithmUserProperties_->value(name);
    }

    ERR("Not a valid property: {:s}", name.toStdString());
    return QVariant();
}

QList<QVariant> VtkAlgorithmProperties::GetUserVectorProperty(QString name) const
{
    if (this->algorithmUserVectorProperties_->contains(name))
    {
        return this->algorithmUserVectorProperties_->value(name);
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
    activeAttributeName_ = name;
}
