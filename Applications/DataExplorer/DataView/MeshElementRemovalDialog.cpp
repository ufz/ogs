/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-16
 * \brief  Implementation of the MeshElementRemovalDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshElementRemovalDialog.h"

#include "Mesh.h"
#include "Elements/Element.h"
#include "MeshLib/Node.h"
#include "MeshEditing/ElementExtraction.h"
#include "AABB.h"
#include "OGSError.h"

#include <QList>
#include <QListWidgetItem>

/// Constructor
MeshElementRemovalDialog::MeshElementRemovalDialog(const ProjectData &project, QDialog* parent)
	: QDialog(parent), _project(project), _currentIndex(0), _aabbIndex(std::numeric_limits<unsigned>::max()), _matIDIndex(std::numeric_limits<unsigned>::max())
{
	setupUi(this);

	const std::vector<MeshLib::Mesh*> mesh_vec (_project.getMeshObjects());

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

MeshElementRemovalDialog::~MeshElementRemovalDialog()
{
}

void MeshElementRemovalDialog::accept()
{
	if (this->newMeshNameEdit->text().size()==0)
	{
		OGSError::box("Please enter name for new mesh.");
		return;
	}

	bool anything_checked (false);

	MeshLib::ElementExtraction ex(*_project.getMesh(this->meshNameComboBox->currentText().toStdString()));
	if (this->elementTypeCheckBox->isChecked())
	{
		QList<QListWidgetItem*> items = this->elementTypeListWidget->selectedItems();
		for (int i=0; i<items.size(); ++i)
			ex.searchByElementType(String2MeshElemType(items[i]->text().toStdString()));
		anything_checked = true;
	}
	if (this->materialIDCheckBox->isChecked())
	{
		QList<QListWidgetItem*> items = this->materialListWidget->selectedItems();
		for (int i=0; i<items.size(); ++i)
			ex.searchByMaterialID(items[i]->text().toInt());
		anything_checked = true;
	}
	if (this->boundingBoxCheckBox->isChecked())
	{
		std::vector<MeshLib::Node*> const& nodes (_project.getMesh(this->meshNameComboBox->currentText().toStdString())->getNodes());
		GeoLib::AABB<MeshLib::Node> const aabb(nodes.begin(), nodes.end());
		MeshLib::Node minAABB = aabb.getMinPoint();
		MeshLib::Node maxAABB = aabb.getMaxPoint();

		// only extract bounding box parameters that have been edited (otherwise there will be rounding errors!)
		minAABB[0] = (aabb_edits[0]) ? this->xMinEdit->text().toDouble() : (minAABB[0]);
		maxAABB[0] = (aabb_edits[1]) ? this->xMaxEdit->text().toDouble() : (maxAABB[0]);
		minAABB[1] = (aabb_edits[2]) ? this->yMinEdit->text().toDouble() : (minAABB[1]);
		maxAABB[1] = (aabb_edits[3]) ? this->yMaxEdit->text().toDouble() : (maxAABB[1]);
		minAABB[2] = (aabb_edits[4]) ? this->zMinEdit->text().toDouble() : (minAABB[2]);
		maxAABB[2] = (aabb_edits[5]) ? this->zMaxEdit->text().toDouble() : (maxAABB[2]);
		ex.searchByBoundingBox(minAABB, maxAABB);
		anything_checked = true;
	}

	if (this->zeroVolumeCheckBox->isChecked())
	{
		ex.searchByContent();
		anything_checked = true;
	}

	if (anything_checked)
	{
		MeshLib::Mesh* new_mesh = ex.removeMeshElements(this->newMeshNameEdit->text().toStdString());
		if (new_mesh)
			emit meshAdded(new_mesh);
		else
		{
			const unsigned error_code (ex.getErrorCode());
			if (error_code == 1)
				OGSError::box("The current selection removes ALL mesh elements.\nPlease change the selection.");
			if (error_code == 2)
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

void MeshElementRemovalDialog::on_boundingBoxCheckBox_toggled(bool is_checked)
{
	this->xMinEdit->setEnabled(is_checked); this->xMaxEdit->setEnabled(is_checked);
	this->yMinEdit->setEnabled(is_checked); this->yMaxEdit->setEnabled(is_checked);
	this->zMinEdit->setEnabled(is_checked); this->zMaxEdit->setEnabled(is_checked);

	if (is_checked && (_currentIndex != _aabbIndex))
	{
		_aabbIndex = _currentIndex;
		std::vector<MeshLib::Node*> const& nodes (_project.getMesh(this->meshNameComboBox->currentText().toStdString())->getNodes());
		GeoLib::AABB<MeshLib::Node> aabb(nodes.begin(), nodes.end());
		MeshLib::Node const minAABB = aabb.getMinPoint();
		MeshLib::Node const maxAABB = aabb.getMaxPoint();
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

void MeshElementRemovalDialog::on_materialIDCheckBox_toggled(bool is_checked)
{
	this->materialListWidget->setEnabled(is_checked);

	if (is_checked && (_currentIndex != _matIDIndex))
	{
		this->materialListWidget->clear();
		_matIDIndex = _currentIndex;
		const std::vector<MeshLib::Element*> elements (_project.getMesh(this->meshNameComboBox->currentText().toStdString())->getElements());
		auto it = std::max_element(elements.begin(), elements.end(),
			[](MeshLib::Element const*const x, MeshLib::Element const*const y) { return x->getValue() < y->getValue(); }
		);
		unsigned max_material ((*it)->getValue());

		for (unsigned i=0; i<=max_material; ++i)
			this->materialListWidget->addItem(QString::number(i));
	}
}

void MeshElementRemovalDialog::on_meshNameComboBox_currentIndexChanged(int idx)
{
	Q_UNUSED(idx);
	this->_currentIndex = this->meshNameComboBox->currentIndex();
	this->newMeshNameEdit->setText(this->meshNameComboBox->currentText() + "_new");
	this->elementTypeListWidget->clearSelection();
	this->materialListWidget->clearSelection();
	if (this->boundingBoxCheckBox->isChecked()) this->on_boundingBoxCheckBox_toggled(true);
	if (this->materialIDCheckBox->isChecked()) this->on_materialIDCheckBox_toggled(true);
}
