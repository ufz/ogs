/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-16
 * \brief  Implementation of the MeshElementRemovalDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshElementRemovalDialog.h"

#include <QList>
#include <QListWidgetItem>
#include <algorithm>

#include "Applications/DataExplorer/Base/OGSError.h"
#include "Applications/DataHolderLib/Project.h"
#include "Elements/Element.h"
#include "GeoLib/AABB.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshLib/Node.h"
#include "MeshLib/Properties.h"
#include "MeshLib/PropertyVector.h"

/// Constructor
MeshElementRemovalDialog::MeshElementRemovalDialog(
    DataHolderLib::Project const& project, QDialog* parent)
    : QDialog(parent),
      _project(project),
      _currentIndex(0),
      _aabbIndex(std::numeric_limits<unsigned>::max()),
      _scalarIndex(std::numeric_limits<unsigned>::max())
{
    setupUi(this);

    auto const& mesh_vec(_project.getMeshObjects());

    const std::size_t nMeshes (mesh_vec.size());
    for (std::size_t i=0; i<nMeshes; ++i)
    {
        std::string name = mesh_vec[i]->getName();
        this->meshNameComboBox->addItem(QString::fromStdString(name));
    }

    if (mesh_vec.empty())
    {
        OGSError::box("No meshes available.");
        QMetaObject::invokeMethod(this, "close", Qt::QueuedConnection);
    }
}

MeshElementRemovalDialog::~MeshElementRemovalDialog() = default;

void MeshElementRemovalDialog::accept()
{
    if (this->newMeshNameEdit->text().size()==0)
    {
        OGSError::box("Please enter name for new mesh.");
        return;
    }

    bool anything_checked (false);

    const MeshLib::Mesh* msh = _project.getMesh(this->meshNameComboBox->currentText().toStdString());
    MeshLib::ElementSearch ex(*msh);
    if (this->elementTypeCheckBox->isChecked())
    {
        QList<QListWidgetItem*> items = this->elementTypeListWidget->selectedItems();
        for (auto& item : items)
            ex.searchByElementType(
                MeshLib::String2MeshElemType(item->text().toStdString()));
        anything_checked = true;
    }
    if (this->scalarArrayCheckBox->isChecked())
    {
        std::string const array_name = this->scalarArrayComboBox->currentText().toStdString();
        double min_val, max_val;
        bool outside = this->outsideButton->isChecked();
        if (outside)
        {
            min_val = this->outsideScalarMinEdit->text().toDouble();
            max_val = this->outsideScalarMaxEdit->text().toDouble();
        }
        else
        {
            min_val = this->insideScalarMinEdit->text().toDouble();
            max_val = this->insideScalarMaxEdit->text().toDouble();
        }
        ex.searchByPropertyValueRange(array_name, min_val, max_val, outside);
        anything_checked = true;
    }
    if (this->boundingBoxCheckBox->isChecked())
    {
        std::vector<MeshLib::Node*> const& nodes (_project.getMesh(this->meshNameComboBox->currentText().toStdString())->getNodes());
        GeoLib::AABB const aabb(nodes.begin(), nodes.end());
        auto minAABB = aabb.getMinPoint();
        auto maxAABB = aabb.getMaxPoint();

        // only extract bounding box parameters that have been edited (otherwise there will be rounding errors!)
        minAABB[0] = (aabb_edits[0]) ? this->xMinEdit->text().toDouble() : (minAABB[0]);
        maxAABB[0] = (aabb_edits[1]) ? this->xMaxEdit->text().toDouble() : (maxAABB[0]);
        minAABB[1] = (aabb_edits[2]) ? this->yMinEdit->text().toDouble() : (minAABB[1]);
        maxAABB[1] = (aabb_edits[3]) ? this->yMaxEdit->text().toDouble() : (maxAABB[1]);
        minAABB[2] = (aabb_edits[4]) ? this->zMinEdit->text().toDouble() : (minAABB[2]);
        maxAABB[2] = (aabb_edits[5]) ? this->zMaxEdit->text().toDouble() : (maxAABB[2]);
        std::vector<MathLib::Point3d> extent;
        extent.push_back(minAABB);
        extent.push_back(maxAABB);
        const GeoLib::AABB updated_aabb(extent.begin(), extent.end());
        ex.searchByBoundingBox(updated_aabb);
        anything_checked = true;
    }

    if (this->zeroVolumeCheckBox->isChecked())
    {
        ex.searchByContent();
        anything_checked = true;
    }

    if (anything_checked)
    {
        MeshLib::Mesh* new_mesh = MeshLib::removeElements(*msh, ex.getSearchedElementIDs(), this->newMeshNameEdit->text().toStdString());
        if (new_mesh)
            emit meshAdded(new_mesh);
        else
        {
            if (new_mesh == nullptr)
                OGSError::box("The current selection removes ALL mesh elements.\nPlease change the selection.");
            if (ex.getSearchedElementIDs().empty())
                OGSError::box("The current selection removes NO mesh elements.");
            delete new_mesh;
            return;
        }
    }
    else
    {
        OGSError::box("No condition set for elements to remove.");
        return;
    }

    this->done(QDialog::Accepted);
}

void MeshElementRemovalDialog::reject()
{
    this->done(QDialog::Rejected);
}

std::size_t MeshElementRemovalDialog::addScalarArrays(MeshLib::Mesh const& mesh) const
{
    MeshLib::Properties const& properties = mesh.getProperties();
    std::vector<std::string> const& names = properties.getPropertyVectorNames();
    for (auto name : names)
    {
        if (properties.existsPropertyVector<int>(name))
        {
            auto const p = properties.getPropertyVector<int>(name);
            if (p->getMeshItemType() == MeshLib::MeshItemType::Cell)
            {
                this->scalarArrayComboBox->addItem(QString::fromStdString(name));
                enableScalarArrayWidgets(true);
            }
        }
        if (properties.existsPropertyVector<double>(name))
        {
            auto const p = properties.getPropertyVector<double>(name);
            if (p->getMeshItemType() == MeshLib::MeshItemType::Cell)
            {
                this->scalarArrayComboBox->addItem(QString::fromStdString(name));
                enableScalarArrayWidgets(true);
            }
        }
    }
    return this->scalarArrayComboBox->count();
}

void MeshElementRemovalDialog::enableScalarArrayWidgets(bool enable) const
{
    this->scalarArrayComboBox->setEnabled(enable);
    this->outsideButton->setEnabled(enable);
    this->insideButton->setEnabled(enable);
    this->outsideScalarMinEdit->setEnabled(enable && this->outsideButton->isChecked());
    this->outsideScalarMaxEdit->setEnabled(enable && this->outsideButton->isChecked());
    this->insideScalarMinEdit->setEnabled(enable && this->insideButton->isChecked());
    this->insideScalarMaxEdit->setEnabled(enable && this->insideButton->isChecked());
}

void MeshElementRemovalDialog::toggleScalarEdits(bool outside) const
{
    this->outsideScalarMinEdit->setEnabled(outside);
    this->outsideScalarMaxEdit->setEnabled(outside);
    this->insideScalarMinEdit->setEnabled(!outside);
    this->insideScalarMaxEdit->setEnabled(!outside);
}

void MeshElementRemovalDialog::on_insideButton_toggled(bool is_checked)
{
    if (this->insideButton->isChecked())
        toggleScalarEdits(false);
    else
        toggleScalarEdits(true);
}

void MeshElementRemovalDialog::on_boundingBoxCheckBox_toggled(bool is_checked)
{
    this->xMinEdit->setEnabled(is_checked); this->xMaxEdit->setEnabled(is_checked);
    this->yMinEdit->setEnabled(is_checked); this->yMaxEdit->setEnabled(is_checked);
    this->zMinEdit->setEnabled(is_checked); this->zMaxEdit->setEnabled(is_checked);

    if (is_checked && (_currentIndex != _aabbIndex))
    {
        _aabbIndex = _currentIndex;
        std::vector<MeshLib::Node*> const& nodes (_project.getMesh(this->meshNameComboBox->currentText().toStdString())->getNodes());
        GeoLib::AABB aabb(nodes.begin(), nodes.end());
        auto const& minAABB = aabb.getMinPoint();
        auto const& maxAABB = aabb.getMaxPoint();
        this->xMinEdit->setText(QString::number(minAABB[0], 'f'));
        this->xMaxEdit->setText(QString::number(maxAABB[0], 'f'));
        this->yMinEdit->setText(QString::number(minAABB[1], 'f'));
        this->yMaxEdit->setText(QString::number(maxAABB[1], 'f'));
        this->zMinEdit->setText(QString::number(minAABB[2], 'f'));
        this->zMaxEdit->setText(QString::number(maxAABB[2], 'f'));
        aabb_edits.fill(false);
    }
}

void MeshElementRemovalDialog::on_elementTypeCheckBox_toggled(bool is_checked)
{
    this->elementTypeListWidget->setEnabled(is_checked);
}

void MeshElementRemovalDialog::on_scalarArrayCheckBox_toggled(bool is_checked)
{
    if (!is_checked)
    {
        enableScalarArrayWidgets(false);
        return;
    }

    MeshLib::Mesh const* const mesh =
        _project.getMesh(meshNameComboBox->currentText().toStdString());
    if (addScalarArrays(*mesh) > 0)
        enableScalarArrayWidgets(true);
    else
    {
        enableScalarArrayWidgets(false);
        OGSError::box("No scalar arrays found");
        on_scalarArrayCheckBox_toggled(false);
    }
}

void MeshElementRemovalDialog::on_meshNameComboBox_currentIndexChanged(int idx)
{
    Q_UNUSED(idx);
    this->_currentIndex = this->meshNameComboBox->currentIndex();
    this->newMeshNameEdit->setText(this->meshNameComboBox->currentText() + "_new");
    this->elementTypeListWidget->clearSelection();
    this->scalarArrayComboBox->clear();
    this->outsideScalarMinEdit->setText("");
    this->outsideScalarMaxEdit->setText("");
    this->insideScalarMinEdit->setText("");
    this->insideScalarMaxEdit->setText("");
    if (this->scalarArrayCheckBox->isChecked())
        on_scalarArrayCheckBox_toggled(true);
    if (this->boundingBoxCheckBox->isChecked())
        on_boundingBoxCheckBox_toggled(true);
}

void MeshElementRemovalDialog::on_scalarArrayComboBox_currentIndexChanged(int idx)
{
    Q_UNUSED(idx);
    MeshLib::Mesh const* const mesh =
        _project.getMesh(meshNameComboBox->currentText().toStdString());
    MeshLib::Properties const& properties = mesh->getProperties();

    std::string const vec_name(scalarArrayComboBox->currentText().toStdString());
    if (vec_name == "")
        return;

    if (properties.existsPropertyVector<int>(vec_name))
    {
        MeshLib::PropertyVector<int> const& vec =
            *properties.getPropertyVector<int>(vec_name);
        auto min = std::min_element(vec.cbegin(), vec.cend());
        auto max = std::max_element(vec.cbegin(), vec.cend());
        this->outsideScalarMinEdit->setText(QString::number(*min));
        this->outsideScalarMaxEdit->setText(QString::number(*max));
        this->insideScalarMinEdit->setText(QString::number(*min));
        this->insideScalarMaxEdit->setText(QString::number(*max));
    }
    else if (properties.existsPropertyVector<double>(vec_name))
    {
        MeshLib::PropertyVector<double> const& vec =
            *properties.getPropertyVector<double>(vec_name);
        auto min = std::min_element(vec.cbegin(), vec.cend());
        auto max = std::max_element(vec.cbegin(), vec.cend());
        this->outsideScalarMinEdit->setText(QString::number(*min));
        this->outsideScalarMaxEdit->setText(QString::number(*max));
        this->insideScalarMinEdit->setText(QString::number(*min));
        this->insideScalarMaxEdit->setText(QString::number(*max));
    }
}
