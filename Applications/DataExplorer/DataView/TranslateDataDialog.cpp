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

    auto const geoNames = geo_models->getGeometryNames();

    std::size_t nGeoObjects(geoNames.size());

    QStringList list;
    for (unsigned i = 0; i < nGeoObjects; ++i)
    {
        this->geoListBox->addItem(QString::fromStdString(geoNames[i]));
    }

    for (int model_index = 0; model_index < mesh_model->rowCount();
         ++model_index)
    {
        auto const* mesh =
            mesh_model->getMesh(mesh_model->index(model_index, 0));
        this->meshListBox->addItem(QString::fromStdString(mesh->getName()));
    }
}

TranslateDataDialog::~TranslateDataDialog() = default;

void TranslateDataDialog::accept()
{
    if (this->meshListBox->count() == 0 and this->geoListBox->count() == 0)
    {
        OGSError::box("Please specify the input data.");
        return;
    }

    QString xinput = this->xlineEdit->text();
    QString yinput = this->ylineEdit->text();
    QString zinput = this->zlineEdit->text();

    if (xinput.isEmpty() || yinput.isEmpty() || zinput.isEmpty())
    {
        OGSError::box("Please specify coordinates to translate the data.");
        return;
    }

    MeshLib::Mesh const* mesh(_mesh_model->getMesh(
        _mesh_model->index(this->meshListBox->currentIndex(), 0)));

    if (mesh == nullptr)
    {
        OGSError::box("No mesh was found.");
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
    if (fabs(xinput.toDouble()) < std::numeric_limits<double>::epsilon() &&
        fabs(yinput.toDouble()) < std::numeric_limits<double>::epsilon() &&
        fabs(zinput.toDouble()) < std::numeric_limits<double>::epsilon())
    {
        GeoLib::AABB aabb(mesh->getNodes().begin(), mesh->getNodes().end());
        auto const [min, max] = aabb.getMinMaxPoints();
        displacement = -(max + min) / 2.0;
    }

    INFO("translate model ({:f}, {:f}, {:f}).",
         displacement[0],
         displacement[1],
         displacement[2]);

    auto const geo_name =
        _geo_models->getGeometryNames()[this->meshListBox->currentIndex()];

    std::vector<GeoLib::Point*> const* point_vec =
        _geo_models->getPointVec(geo_name);
    std::size_t const n_points = point_vec->size();
    for (std::size_t i = 0; i < n_points; ++i)
    {
        for (std::size_t c = 0; c < 3; ++c)
        {
            (*(*point_vec)[i])[c] += displacement[c];
        }
    }

    MeshLib::moveMeshNodes(mesh->getNodes().begin(), mesh->getNodes().end(),
                           displacement);

    _mesh_model->updateMesh(const_cast<MeshLib::Mesh*>(mesh));

    _geo_models->updateGeometry(geo_name);

    this->done(QDialog::Accepted);
}
