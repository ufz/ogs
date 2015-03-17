/**
 * \file
 * \author Lars Bilke
 * \date   2009-09-24
 * \brief  Definition of the MshView class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#ifndef MSHVIEW_H
#define MSHVIEW_H

#include "Point.h"
#include "GeoType.h"
#include <QTreeView>

class MshModel;
class VtkMeshSource;
class vtkUnstructuredGridAlgorithm;

namespace MeshLib {
	class Mesh;
}
/*
   namespace GeoLib {
    class Point;
   }
 */
/**
 *	The DataView is table view which acts as a base class for displaying
 *  several OSG data formats.
 */
class MshView : public QTreeView
{
	Q_OBJECT

public:
	MshView(QWidget* parent = 0);
	~MshView();

public slots:
	void updateView();

protected slots:
	/// Is called when the selection of this view changes.
	void selectionChanged(const QItemSelection &selected, const QItemSelection &deselected);

private:
	void contextMenuEvent( QContextMenuEvent* event );

private slots:
	/// Opens a dialog for editing meshes.
	void openMeshEditDialog();

	/// Opens a dialog for editing material groups.
	void openValuesEditDialog();

	/// Adds a new mesh.
	void addMesh();

	void addDIRECTSourceTerms();

	void extractSurfaceMesh();

	void exportToTetGen();

	void loadDIRECTSourceTerms();

	void convertMeshToGeometry();

	void exportToShapefile() const;

	/// Remove the currently selected mesh.
	void removeMesh();

	/// Calls the FileDialog to save a mesh to a file.
	void writeToFile() const;

	/// Calls the dialog for calculating an element quality metric
	void checkMeshQuality();

signals:
	void elementSelected(vtkUnstructuredGridAlgorithm const*const, unsigned, bool);
	void enableSaveButton(bool);
	void enableRemoveButton(bool);
	void meshSelected(MeshLib::Mesh const*const);
	void openMeshFile(int);
	void qualityCheckRequested(vtkUnstructuredGridAlgorithm*);
	void removeSelectedMeshComponent();
	void requestCondSetupDialog(const std::string&, const GeoLib::GEOTYPE, const std::size_t, bool on_points);
	void requestMeshRemoval(const QModelIndex&);
	void requestMeshToGeometryConversion(const MeshLib::Mesh*);
	void loadFEMCondFileRequested(const std::string);
	void saveMeshAction();

};
#endif // MSHVIEW_H
