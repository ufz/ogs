/**
 * \file
 * \author Karsten Rink
 * \date   2010-12-09
 * \brief  Implementation of the LineEditDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LineEditDialog.h"
#include "OGSError.h"
#include <QStringList>
#include <QStringListModel>

LineEditDialog::LineEditDialog(const GeoLib::PolylineVec &ply_vec, QDialog* parent)
    : QDialog(parent), _allPly(new QStringListModel), _selPly(new QStringListModel),
      _geoName(ply_vec.getName())
{
    setupUi(this);

    this->proximityEdit->setValidator(new QDoubleValidator(0, 100, 8, this));

    std::size_t nPly(ply_vec.size());
    QStringList list;
    for (std::size_t i = 0; i < nPly; i++)
    {
        std::string ply_name("");
        ply_vec.getNameOfElementByID(i, ply_name);
        list.append("Line " + QString::number(i) + "  " + QString::fromStdString(ply_name));
    }
    _allPly->setStringList(list);

    this->allPlyView->setModel(_allPly);
    this->selectedPlyView->setModel(_selPly);
}

LineEditDialog::~LineEditDialog()
{
    delete _allPly;
    delete _selPly;
}

void LineEditDialog::on_selectPlyButton_pressed()
{
    QModelIndexList selected = this->allPlyView->selectionModel()->selectedIndexes();
    QStringList list = _selPly->stringList();

    for (auto& it : selected)
    {
        list.append(it.data().toString());

        _allPly->removeRow(it.row());
    }
    _selPly->setStringList(list);
}

void LineEditDialog::on_deselectPlyButton_pressed()
{
    QModelIndexList selected = this->selectedPlyView->selectionModel()->selectedIndexes();
    QStringList list = _allPly->stringList();

    for (auto& it : selected)
    {
        list.append(it.data().toString());

        _selPly->removeRow(it.row());
    }
    _allPly->setStringList(list);
}

void LineEditDialog::accept()
{
    std::vector<std::size_t> selectedIndeces = this->getSelectedIndeces(_selPly->stringList());

    if (!selectedIndeces.empty())
    {
        std::string prox_string = this->proximityEdit->text().toStdString();
        double prox =
            (prox_string.empty()) ? 0.0 : strtod(prox_string.c_str(), nullptr);
        std::string ply_name =
                (plyNameEdit->text().toStdString().empty()) ? "" : plyNameEdit->text().
                toStdString();
        emit connectPolylines(_geoName,
                              selectedIndeces,
                              prox,
                              ply_name,
                              this->closePlyCheckBox->isChecked(),
                              this->createSfcCheckBox->isChecked());
        this->done(QDialog::Accepted);
    }
    else
        OGSError::box("No polylines selected", "Error");
}

void LineEditDialog::reject()
{
    this->done(QDialog::Rejected);
}

std::vector<std::size_t> LineEditDialog::getSelectedIndeces(QStringList list)
{
    std::vector<std::size_t> indexList;
    for (auto& it : list)
    {
        QString s = it.mid(5, it.indexOf("  ") - 5);
        indexList.push_back(s.toInt());
    }
    return indexList;
}
