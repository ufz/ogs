/**
 * \file
 * \author Karsten Rink
 * \date   2010-12-09
 * \brief  Implementation of the LineEditDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    : QDialog(parent), allPly_(new QStringListModel), selPly_(new QStringListModel),
      geoName_(ply_vec.getName())
{
    setupUi(this);

    this->proximityEdit->setValidator(new QDoubleValidator(0, 100, 8, this));

    std::size_t nPly(ply_vec.size());
    QStringList list;
    for (std::size_t i = 0; i < nPly; i++)
    {
        std::string ply_name;
        ply_vec.getNameOfElementByID(i, ply_name);
        list.append("Line " + QString::number(i) + "  " + QString::fromStdString(ply_name));
    }
    allPly_->setStringList(list);

    this->allPlyView->setModel(allPly_);
    this->selectedPlyView->setModel(selPly_);
}

LineEditDialog::~LineEditDialog()
{
    delete allPly_;
    delete selPly_;
}

void LineEditDialog::on_selectPlyButton_pressed()
{
    QModelIndexList selected = this->allPlyView->selectionModel()->selectedIndexes();
    QStringList list = selPly_->stringList();

    for (auto& index : selected)
    {
        list.append(index.data().toString());

        allPly_->removeRow(index.row());
    }
    selPly_->setStringList(list);
}

void LineEditDialog::on_deselectPlyButton_pressed()
{
    QModelIndexList selected = this->selectedPlyView->selectionModel()->selectedIndexes();
    QStringList list = allPly_->stringList();

    for (auto& index : selected)
    {
        list.append(index.data().toString());

        selPly_->removeRow(index.row());
    }
    allPly_->setStringList(list);
}

void LineEditDialog::accept()
{
    std::vector<std::size_t> selectedIndeces = this->getSelectedIndeces(selPly_->stringList());

    if (!selectedIndeces.empty())
    {
        std::string prox_string = this->proximityEdit->text().toStdString();
        double prox =
            (prox_string.empty()) ? 0.0 : strtod(prox_string.c_str(), nullptr);
        std::string ply_name =
                (plyNameEdit->text().toStdString().empty()) ? "" : plyNameEdit->text().
                toStdString();
        emit connectPolylines(geoName_,
                              selectedIndeces,
                              prox,
                              ply_name,
                              this->closePlyCheckBox->isChecked(),
                              this->createSfcCheckBox->isChecked());
        this->done(QDialog::Accepted);
    }
    else
    {
        OGSError::box("No polylines selected", "Error");
    }
}

void LineEditDialog::reject()
{
    this->done(QDialog::Rejected);
}

std::vector<std::size_t> LineEditDialog::getSelectedIndeces(QStringList list)
{
    std::vector<std::size_t> indexList;
    for (auto& index : list)
    {
        QString s = index.mid(5, index.indexOf("  ") - 5);
        indexList.push_back(s.toInt());
    }
    return indexList;
}
