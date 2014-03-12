/**
 * \file   MeshAnalysisDialog.h
 * \author Karsten Rink
 * \date   2014-02-24
 * \brief  Definition of the MeshAnalysisDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHANALYSISDIALOG_H
#define MESHANALYSISDIALOG_H

#include "ui_MeshAnalysis.h"
#include <QDialog>

#include "MeshQuality/ElementErrorCode.h"

namespace MeshLib {
	class Mesh;
}

/**
 * \brief A dialog window for calling mesh analysis methods
 */
class MeshAnalysisDialog : public QDialog, private Ui_MeshAnalysis
{
	Q_OBJECT

public:
	MeshAnalysisDialog(const std::vector<MeshLib::Mesh*> &mesh_vec, QDialog* parent = nullptr);
	~MeshAnalysisDialog(void);

private:
	/// Prepares the output for the node message window
	void nodesMsgOutput(const std::vector<std::size_t> &node_ids, unsigned nCollapsableNodes);

	const std::vector<MeshLib::Mesh*>& _mesh_vec;

private slots:
	/// Starts the analysis
	void on_startButton_pressed();

	/// Closes the dialog
	void on_closeButton_pressed() { this->close(); }
};

#endif //MESHANALYSISDIALOG_H
