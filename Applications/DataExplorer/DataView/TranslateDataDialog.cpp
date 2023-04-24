/**
 * \file
 * \author Karsten Rink
 * \date   2014-02-24
 * \brief  Implementation of the TranslateDataDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TranslateDataDialog.h"

#include <QStringList>

#include "Base/StrictDoubleValidator.h"
#include "GEOModels.h"
#include "GeoLib/AABB.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/moveMeshNodes.h"
#include "MeshLib/Node.h"
#include "MeshModel.h"

TranslateDataDialog::TranslateDataDialog(MeshModel* mesh_model,
                                         GEOModels* geo_models,
                                         QDialog* parent)
    : QDialog(parent), _mesh_model(mesh_model), _geo_models(geo_models)
{
    setupUi(this);
    assert(_geo_models != nullptr);
    assert(_mesh_model != nullptr);
    auto const geoNames = _geo_models->getGeometryNames();

    QStringList dataList;
    for (auto const& name : geoNames)
    {
        dataList.append(QString::fromStdString(name));
    }

    for (int model_index = 0; model_index < _mesh_model->rowCount();
         ++model_index)
    {
        auto const* mesh =
            _mesh_model->getMesh(_mesh_model->index(model_index, 0));
        dataList.append(QString::fromStdString(mesh->getName()));
    }

    if (dataList.empty())
    {
        this->selectDataButton->setDisabled(true);
        this->deselectDataButton->setDisabled(true);
        dataList.append("[No data available.]");
    }

    _allData.setStringList(dataList);
    this->allDataView->setModel(&_allData);
    this->selectedDataView->setModel(&_selData);
}

void TranslateDataDialog::on_selectDataButton_pressed()
{
    QModelIndexList const selected =
        this->allDataView->selectionModel()->selectedIndexes();
    QStringList list = _selData.stringList();

    for (auto& index : selected)
    {
        list.append(index.data().toString());

        _allData.removeRow(index.row());
    }
    _selData.setStringList(list);
}

void TranslateDataDialog::on_deselectDataButton_pressed()
{
    QModelIndexList selected =
        this->selectedDataView->selectionModel()->selectedIndexes();
    QStringList list = _allData.stringList();

    for (auto& index : selected)
    {
        list.append(index.data().toString());

        _selData.removeRow(index.row());
    }
    _allData.setStringList(list);
}

void TranslateDataDialog::moveGeometry(Eigen::Vector3d displacement,
                                       const std::string name)
{
    std::vector<GeoLib::Point*> const* point_vec =
        _geo_models->getPointVec(name);
    if (point_vec == nullptr)
    {
        OGSError::box("The geometry is faulty.");
        return;
    }
    std::size_t const n_points = point_vec->size();
    for (std::size_t i = 0; i < n_points; ++i)
    {
        for (std::size_t c = 0; c < 3; ++c)
        {
            (*(*point_vec)[i])[c] += displacement[c];
        }
    }
    for (auto* point : *point_vec)
    {
        point->asEigenVector3d() += displacement;
    }
    _geo_models->updateGeometry(name);
}

void TranslateDataDialog::moveMesh(Eigen::Vector3d displacement,
                                   const std::string name)
{
    MeshLib::Mesh const* mesh(_mesh_model->getMesh(name));
    if (mesh == nullptr)
    {
        OGSError::box("The mesh is faulty.");
        return;
    }
    /*
    if (std::fabs(xinput.toDouble()) < std::numeric_limits<double>::epsilon() &&
        std::fabs(yinput.toDouble()) < std::numeric_limits<double>::epsilon() &&
        std::fabs(zinput.toDouble()) < std::numeric_limits<double>::epsilon())
    {
        GeoLib::AABB const aabb(mesh->getNodes().begin(),
                                mesh->getNodes().end());
        auto const [min, max] = aabb.getMinMaxPoints();
        displacement = -(max + min) / 2.0;
    }
    */
    MeshLib::moveMeshNodes(mesh->getNodes().begin(), mesh->getNodes().end(),
                           displacement);
    _mesh_model->updateMesh(const_cast<MeshLib::Mesh*>(mesh));
}

void TranslateDataDialog::accept()
{
    if (this->_selData.rowCount() == 0)
    {
        OGSError::box("Please specify the input data.");
        return;
    }

    QString const xinput = this->xlineEdit->text();
    QString const yinput = this->ylineEdit->text();
    QString const zinput = this->zlineEdit->text();

    if (xinput.isEmpty() || yinput.isEmpty() || zinput.isEmpty())
    {
        OGSError::box("Please specify coordinates to translate the data.");
        return;
    }

    // why is 0 not a valid input ?
    if (!xinput.toDouble() or !yinput.toDouble() or !zinput.toDouble())
    {
        OGSError::box("The input must be a real number.");
        return;
    }

    Eigen::Vector3d displacement{xinput.toDouble(), yinput.toDouble(),
                                 zinput.toDouble()};

    INFO("translate model ({:f}, {:f}, {:f}).",
         displacement[0],
         displacement[1],
         displacement[2]);

    std::vector<std::string> const selectedData =
        this->getSelectedObjects(_selData.stringList());

    auto const geoNames = _geo_models->getGeometryNames();

    for (auto const& data_name: selectedData)
    {
        if (std::find(std::begin(geoNames), std::end(geoNames), data_name) !=
            std::end(geoNames))
        {
            moveGeometry(displacement, data_name);
            continue;
        }
        moveMesh(displacement, data_name);
    }

    this->done(QDialog::Accepted);
}

std::vector<std::string> TranslateDataDialog::getSelectedObjects(
    QStringList list)
{
    std::vector<std::string> indexList;
    std::transform(list.begin(), list.end(), std::back_inserter(indexList),
                   [](auto const& index) { return index.toStdString(); });
    return indexList;
}