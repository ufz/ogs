/**
 * \file   MeshAnalysisDialog.cpp
 * \author Karsten Rink
 * \date   2014-02-24
 * \brief  Implementation of the MeshAnalysisDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshAnalysisDialog.h"
#include "Mesh.h"
#include "MeshQuality/MeshValidation.h"
#include "MeshEditing/MeshRevision.h"
#include "MeshSearch/NodeSearch.h"

#include "StrictDoubleValidator.h"

#include <logog/include/logog.hpp>

MeshAnalysisDialog::MeshAnalysisDialog(
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& mesh_vec,
    QDialog* parent)
: QDialog(parent), _mesh_vec(mesh_vec)
{
    setupUi(this);

    if (mesh_vec.empty())
        this->startButton->setDisabled(true);

    for (std::size_t i=0; i<mesh_vec.size(); ++i)
        this->meshListBox->addItem(QString::fromStdString(mesh_vec[i]->getName()));

    StrictDoubleValidator* collapse_threshold_validator = new StrictDoubleValidator(0, 1000000, 7, this);
    this->collapsibleNodesThreshold->setValidator (collapse_threshold_validator);

    StrictDoubleValidator* volume_threshold_validator = new StrictDoubleValidator(0, 1e10, 10, this);
    this->zeroVolumeThreshold->setValidator (volume_threshold_validator);
}

MeshAnalysisDialog::~MeshAnalysisDialog()
{
}

void MeshAnalysisDialog::on_startButton_pressed()
{
    MeshLib::Mesh const& mesh (*_mesh_vec[this->meshListBox->currentIndex()].get());

    MeshLib::NodeSearch ns(mesh);
    ns.searchUnused();
    const std::vector<std::size_t> unusedNodesIdx (ns.getSearchedNodeIDs());
    MeshLib::MeshRevision rev(const_cast<MeshLib::Mesh&>(mesh));
    std::vector<std::size_t> const& collapsibleNodeIds (rev.collapseNodeIndices(
        this->collapsibleNodesThreshold->text().toDouble() + std::numeric_limits<double>::epsilon()));
    this->nodesGroupBox->setTitle("Nodes (out of " + QString::number(mesh.getNumberOfNodes()) + ")");
    this->nodesMsgOutput(unusedNodesIdx, collapsibleNodeIds);

    const std::vector<ElementErrorCode> element_error_codes (MeshLib::MeshValidation::testElementGeometry(
        mesh, this->zeroVolumeThreshold->text().toDouble() + std::numeric_limits<double>::epsilon()));
    this->elementsGroupBox->setTitle("Elements (out of " + QString::number(mesh.getNumberOfElements()) + ")");
    this->elementsMsgOutput(element_error_codes);

    unsigned const n_holes (MeshLib::MeshValidation::detectHoles(mesh));
    if (n_holes>0)
        this->meshHoleOutputLabel->setText("<strong>" + QString::number(n_holes) + " hole(s) found within the mesh</strong>");
}

void MeshAnalysisDialog::nodesMsgOutput(std::vector<std::size_t> const& node_ids, std::vector<std::size_t> const& collapsibleNodeIds)
{
    const std::size_t nNodeIds (node_ids.size());
    QString nodes_output("");
    if (node_ids.empty())
        nodes_output += "No unused nodes found.";
    else
    {
        nodes_output += (QString::number(nNodeIds) + " nodes are not part of any element:\n");
        for (std::size_t i=0; i<nNodeIds; ++i)
            nodes_output += (QString::number(node_ids[i]) + ", ");
    }
    this->unusedNodesText->setText(nodes_output);

    std::size_t const nNodes(collapsibleNodeIds.size());
    QString node_ids_str("");
    unsigned count(0);
    for (std::size_t i = 0; i < nNodes; ++i)
        if (i != collapsibleNodeIds[i])
        {
            node_ids_str.append(QString::number(i) + ", ");
            count++;
        }
    nodes_output = (count > 0) ? QString::number(count) + " nodes found:\n" : "No nodes found.";
    nodes_output.append(node_ids_str);
    this->collapsibleNodesText->setText(nodes_output);
}

void MeshAnalysisDialog::elementsMsgOutput(const std::vector<ElementErrorCode> &element_error_codes)
{
    std::array<std::string, static_cast<std::size_t>(ElementErrorFlag::MaxValue)> output_str(MeshLib::MeshValidation::ElementErrorCodeOutput(element_error_codes));

    this->zeroVolumeText->setText(QString::fromStdString(output_str[0]));
    this-> nonPlanarText->setText(QString::fromStdString(output_str[1]));
    this-> nonConvexText->setText(QString::fromStdString(output_str[2]));
    this-> nodeOrderText->setText(QString::fromStdString(output_str[3]));
}
