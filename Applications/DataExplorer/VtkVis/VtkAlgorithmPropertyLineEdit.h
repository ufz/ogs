/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-18
 * \brief  Definition of the VtkAlgorithmPropertyLineEdit class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QLineEdit>
#include <QVariant>

class VtkAlgorithmProperties;
class QString;

/// @brief This QLineEdit sets a user property on the given VtkAlgorithmProperties
/// object automatically.
class VtkAlgorithmPropertyLineEdit : public QLineEdit
{
    Q_OBJECT

public:
    /// @brief Constructor.
    /// @param contents The initial text.
    /// @param name The name of the user property to set.
    /// @param type The type of the property.
    /// @param algProps The VtkAlgorithmProperties object.
    /// @param parent The parent widget.
    VtkAlgorithmPropertyLineEdit(const QString& contents,
                                 QString name,
                                 QVariant::Type type,
                                 VtkAlgorithmProperties* algProps,
                                 QWidget* parent = nullptr);
    ~VtkAlgorithmPropertyLineEdit() override;

private:
    const QString _name;
    VtkAlgorithmProperties* _algProps;
    QVariant::Type _type;

private slots:
    /// @brief This slots is automatically called when the text changed.
    void setNewValue();
};
