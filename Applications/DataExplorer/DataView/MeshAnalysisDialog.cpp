/**
 * \file
 * \author Karsten Rink
 * \date   2014-02-24
 * \brief  Implementation of the MeshAnalysisDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshAnalysisDialog.h"

#include "Base/StrictDoubleValidator.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshToolsLib/MeshEditing/MeshRevision.h"
#include "MeshToolsLib/MeshQuality/MeshValidation.h"

MeshAnalysisDialog::MeshAnalysisDialog(
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& mesh_vec,
    QDialog* parent)
    : QDialog(parent), _mesh_vec(mesh_vec)
{
    setupUi(this);

    if (mesh_vec.empty())
    {
        this->startButton->setDisabled(true);
    }

    for (const auto& mesh : mesh_vec)
    {
        this->meshListBox->addItem(QString::fromStdString(mesh->getName()));
    }

    auto* collapse_threshold_validator =
        new StrictDoubleValidator(0, 1000000, 7, this);
    this->collapsibleNodesThreshold->setValidator(collapse_threshold_validator);

    auto* volume_threshold_validator =
        new StrictDoubleValidator(0, 1e10, 10, this);
    this->zeroVolumeThreshold->setValidator(volume_threshold_validator);
}

MeshAnalysisDialog::~MeshAnalysisDialog() = default;

void MeshAnalysisDialog::on_startButton_pressed()
{
    MeshLib::Mesh const& mesh(
        *_mesh_vec[this->meshListBox->currentIndex()].get());

    MeshLib::NodeSearch ns(mesh);
    ns.searchUnused();
    const std::vector<std::size_t> unusedNodesIdx(ns.getSearchedNodeIDs());
    MeshToolsLib::MeshRevision rev(const_cast<MeshLib::Mesh&>(mesh));
    std::vector<std::size_t> const& collapsibleNodeIds(rev.collapseNodeIndices(
        this->collapsibleNodesThreshold->text().toDouble() +
        std::numeric_limits<double>::epsilon()));
    this->nodesGroupBox->setTitle(
        "Nodes (out of " + QString::number(mesh.getNumberOfNodes()) + ")");
    this->nodesMsgOutput(unusedNodesIdx, collapsibleNodeIds);

    const std::vector<ElementErrorCode> element_error_codes(
        MeshToolsLib::MeshValidation::testElementGeometry(
            mesh,
            this->zeroVolumeThreshold->text().toDouble() +
                std::numeric_limits<double>::epsilon()));
    this->elementsGroupBox->setTitle(
        "Elements (out of " + QString::number(mesh.getNumberOfElements()) +
        ")");
    this->elementsMsgOutput(element_error_codes);

    unsigned const n_holes(MeshToolsLib::MeshValidation::detectHoles(mesh));
    if (n_holes > 0)
    {
        this->meshHoleOutputLabel->setText(
            "<strong>" + QString::number(n_holes) +
            " hole(s) found within the mesh</strong>");
    }
}

void MeshAnalysisDialog::nodesMsgOutput(
    std::vector<std::size_t> const& node_ids,
    std::vector<std::size_t> const& collapsibleNodeIds)
{
    const std::size_t nNodeIds(node_ids.size());
    QString nodes_output("");
    if (node_ids.empty())
    {
        nodes_output += "No unused nodes found.";
    }
    else
    {
        nodes_output += (QString::number(nNodeIds) +
                         " nodes are not part of any element:\n");
        for (std::size_t i = 0; i < nNodeIds; ++i)
        {
            nodes_output += (QString::number(node_ids[i]) + ", ");
        }
    }
    this->unusedNodesText->setText(nodes_output);

    std::size_t const nNodes(collapsibleNodeIds.size());
    QString node_ids_str("");
    unsigned count(0);
    for (std::size_t i = 0; i < nNodes; ++i)
    {
        if (i != collapsibleNodeIds[i])
        {
            node_ids_str.append(QString::number(i) + ", ");
            count++;
        }
    }
    nodes_output = (count > 0) ? QString::number(count) + " nodes found:\n"
                               : "No nodes found.";
    nodes_output.append(node_ids_str);
    this->collapsibleNodesText->setText(nodes_output);
}

void MeshAnalysisDialog::elementsMsgOutput(
    const std::vector<ElementErrorCode>& element_error_codes)
{
    std::array<std::string,
               static_cast<std::size_t>(ElementErrorFlag::MaxValue)>
        output_str(MeshToolsLib::MeshValidation::ElementErrorCodeOutput(
            element_error_codes));

    this->zeroVolumeText->setText(QString::fromStdString(output_str[0]));
    this->nonPlanarText->setText(QString::fromStdString(output_str[1]));
    this->nonConvexText->setText(QString::fromStdString(output_str[2]));
    this->nodeOrderText->setText(QString::fromStdString(output_str[3]));
}
