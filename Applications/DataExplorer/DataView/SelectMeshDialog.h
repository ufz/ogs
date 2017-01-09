/**
 * \file
 * \author Karsten Rink
 * \date   2012-04-20
 * \brief  Definition of the SelectMeshDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SELECTMESHDIALOG_H
#define SELECTMESHDIALOG_H

#include <QDialog>

namespace GeoLib {
    class GeoObject;
}

class QDialogButtonBox;
class QLabel;
class QComboBox;
class QVBoxLayout;

/**
 * \brief Small dialog for setting a name for an object.
 */
class SelectMeshDialog : public QDialog
{
    Q_OBJECT

public:
    /// Constructor
    SelectMeshDialog(const GeoLib::GeoObject* geo_object,
                  const std::list<std::string> &msh_names,
                  QDialog* parent = 0);
    ~SelectMeshDialog();

    QDialogButtonBox* _buttonBox; /// The buttons used in this dialog.

private:
    /// Constructs a dialog window
    void setupDialog(const std::list<std::string> &msh_names);

    QLabel* _txt_label;
    QComboBox* _msh_names;
    QVBoxLayout* _layout;
    const GeoLib::GeoObject* _geo_object;


private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept();

    /// Instructions if the Cancel-Button has been pressed.
    void reject();

signals:
    //void requestNameChange(const std::string&, const GeoLib::GEOTYPE, std::size_t, std::string);
};

#endif //SELECTMESHDIALOG_H
