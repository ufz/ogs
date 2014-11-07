/**
 * \file   MeshAnalysisDialog
 * \author Karsten Rink
 * \date   2014-02-24
 * \brief  Implementation of the MeshAnalysisDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshAnalysisDialog.h"
#include "Mesh.h"
#include "MeshQuality/MeshValidation.h"
#include "MeshEditing/MeshRevision.h"

#include "StrictDoubleValidator.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"


MeshAnalysisDialog::MeshAnalysisDialog(const std::vector<MeshLib::Mesh*> &mesh_vec, QDialog* parent)
	: QDialog(parent), _mesh_vec (mesh_vec)
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
	const MeshLib::Mesh* mesh (_mesh_vec[this->meshListBox->currentIndex()]);

	const std::vector<std::size_t> unusedNodesIdx (MeshLib::MeshValidation::removeUnusedMeshNodes(*const_cast<MeshLib::Mesh*>(mesh)));
	MeshLib::MeshRevision rev(const_cast<MeshLib::Mesh&>(*mesh));
	std::vector<std::size_t> const& collapsibleNodeIds (rev.collapseNodeIndeces(
		this->collapsibleNodesThreshold->text().toDouble() + std::numeric_limits<double>::epsilon()));
	this->nodesGroupBox->setTitle("Nodes (out of " + QString::number(mesh->getNNodes()) + ")");
	this->nodesMsgOutput(unusedNodesIdx, collapsibleNodeIds);

	const std::vector<ElementErrorCode> element_error_codes (MeshLib::MeshValidation::testElementGeometry(
		*mesh, this->zeroVolumeThreshold->text().toDouble() + std::numeric_limits<double>::epsilon()));
	this->elementsGroupBox->setTitle("Elements (out of " + QString::number(mesh->getNElements()) + ")");
	this->elementsMsgOutput(element_error_codes);
}

void MeshAnalysisDialog::nodesMsgOutput(std::vector<std::size_t> const& node_ids, std::vector<std::size_t> const& collapsibleNodeIds)
{
	const std::size_t nNodeIds (node_ids.size());
	QString nodes_output("");
	if (node_ids.empty())
		nodes_output += "No unused nodes found.";
	else
	{
		(QString::number(nNodeIds) + " nodes are not part of any element:\n");
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
