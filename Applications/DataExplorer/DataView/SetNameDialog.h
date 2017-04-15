/**
 * \file
 * \author Karsten Rink
 * \date   2011-10-26
 * \brief  Definition of the SetNameDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "GeoType.h"

#include <QDialog>

class QDialogButtonBox;
class QLabel;
class QLineEdit;
class QVBoxLayout;

/**
 * \brief Small dialog for setting a name for an object.
 */
class SetNameDialog : public QDialog
{
    Q_OBJECT

public:
    /// Constructor
    SetNameDialog(const std::string& geo_object_type, std::size_t id,
                  const std::string& old_name, QDialog* parent = nullptr);
    ~SetNameDialog();

    std::string getNewName();

private:
    /// Constructs a dialog window
    void setupDialog(const QString &label, const std::string &old_name);

    QLabel* _txt_label;
    QLineEdit* _new_name;
    QVBoxLayout* _layout;
    QDialogButtonBox* _buttonBox;

private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept();

    /// Instructions if the Cancel-Button has been pressed.
    void reject();
};
