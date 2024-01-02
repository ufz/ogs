/**
 * \file
 * \author Karsten Rink
 * \date   2014-02-24
 * \brief  Implementation of the TranslateDataDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TranslateDataDialog.h"

#include <QStringList>

#include "Base/StrictDoubleValidator.h"
#include "Base/Utils.h"
#include "GEOModels.h"
#include "GeoLib/AABB.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshModel.h"
#include "MeshToolsLib/MeshEditing/moveMeshNodes.h"

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
    Utils::moveSelectedItems(this->allDataView, _allData, _selData);
}

void TranslateDataDialog::on_deselectDataButton_pressed()
{
    Utils::moveSelectedItems(this->selectedDataView, _selData, _allData);
}

void TranslateDataDialog::moveGeometry(Eigen::Vector3d const& displacement,
                                       std::string const& name)
{
    std::vector<GeoLib::Point*> const* point_vec =
        _geo_models->getPointVec(name);
    if (point_vec == nullptr)
    {
        OGSError::box("The geometry is faulty.");
        return;
    }

    MeshToolsLib::moveMeshNodes(point_vec->begin(), point_vec->end(),
                                displacement);

    _geo_models->updateGeometry(name);
}

void TranslateDataDialog::moveMesh(Eigen::Vector3d const& displacement,
                                   std::string const& name)
{
    MeshLib::Mesh const* mesh(_mesh_model->getMesh(name));
    if (mesh == nullptr)
    {
        OGSError::box("The mesh is faulty.");
        return;
    }

    MeshToolsLib::moveMeshNodes(mesh->getNodes().begin(),
                                mesh->getNodes().end(), displacement);
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

    bool ok;
    if (!xinput.toDouble(&ok) or !yinput.toDouble(&ok) or !zinput.toDouble(&ok))
    {
        INFO(
            "If the x/y/z-input is 0, not specified or not a real number, it "
            "is used as 0.");
    }

    Eigen::Vector3d const displacement{xinput.toDouble(), yinput.toDouble(),
                                       zinput.toDouble()};

    INFO("translate model ({:f}, {:f}, {:f}).",
         displacement[0],
         displacement[1],
         displacement[2]);

    std::vector<std::string> const selectedData =
        Utils::getSelectedObjects(_selData.stringList());

    auto const geoNames = _geo_models->getGeometryNames();

    for (auto const& data_name : selectedData)
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