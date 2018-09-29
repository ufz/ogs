/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-20
 * \brief  Definition of the VtkAlgorithmPropertyCheckbox class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QCheckBox>

class VtkAlgorithmProperties;

/// @brief This checkbox sets a user property on the given VtkAlgorithmProperties
/// object automatically.
class VtkAlgorithmPropertyCheckbox : public QCheckBox
{
    Q_OBJECT

public:
    /// @brief Constructor.
    /// @param value The initial check state.
    /// @param name The name of the user property to set.
    /// @param algProps The VtkAlgorithmProperties object.
    /// @param parent The parent widget.
    VtkAlgorithmPropertyCheckbox(const bool value, QString name,
                                 VtkAlgorithmProperties* algProps,
                                 QWidget* parent = nullptr);

    /// @brief Destructor.
    ~VtkAlgorithmPropertyCheckbox() override;

private:
    const QString _name;
    VtkAlgorithmProperties* _algProps;

private slots:
    /// @brief This slots is automatically called when the checkbox state changed.
    void setNewValue(int state);
};
