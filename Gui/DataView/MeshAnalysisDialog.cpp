/**
 * \file   MeshAnalysisDialog
 * \author Karsten Rink
 * \date   2014-02-24
 * \brief  Implementation of the MeshAnalysisDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshAnalysisDialog.h"
#include "Mesh.h"
#include "MeshQuality/MeshValidation.h"

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
}

MeshAnalysisDialog::~MeshAnalysisDialog()
{
}

void MeshAnalysisDialog::on_startButton_pressed()
{
	const MeshLib::Mesh* mesh (_mesh_vec[this->meshListBox->currentIndex()]);

	const std::vector<std::size_t> unusedNodesIdx (MeshLib::MeshValidation::removeUnusedMeshNodes(*const_cast<MeshLib::Mesh*>(mesh)));
	this->nodesMsgOutput(unusedNodesIdx);

	const std::vector<ElementErrorCode> element_error_codes (MeshLib::MeshValidation::testElementGeometry(*mesh));
	//this->elementMsgOutput(element_error_codes);
	QString elementMsgOutput = QString::fromStdString(MeshLib::MeshValidation::ElementErrorCodeOutput(element_error_codes));
	this->elementsMsg->setText(elementMsgOutput);
}

void MeshAnalysisDialog::nodesMsgOutput(const std::vector<std::size_t> &node_ids)
{
	const std::size_t nNodeIds (node_ids.size());
	QString nodes_output("");
	if (node_ids.empty())
		nodes_output += "Nothing to report.";
	else
	{
		(QString::number(nNodeIds) + " nodes are not part of any element:\n");
		for (std::size_t i=0; i<nNodeIds; ++i)
			nodes_output += (QString::number(node_ids[i]) + ", ");
	}
	this->nodesMsg->setText(nodes_output);
}

