/**
 * \file
 * \author Lars Bilke
 * \date   2009-09-24
 * \brief  Definition of the MshView class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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

	/// Selects items without sending signals.
	//void selectionChangedFromOutside(const QItemSelection &selected,
	//	const QItemSelection &deselected);

	/// Clears the selection
	//void clearSelection();

private:
	void contextMenuEvent( QContextMenuEvent* event );

private slots:
	/// Opens a dialog for editing meshes.
	void openMshEditDialog();

	/// Opens a dialog for editing material groups.
	void openValuesEditDialog();

	/// Adds a new mesh.
	void addMesh();

	void addDIRECTSourceTerms();

	void extractSurfaceMesh();

	void loadDIRECTSourceTerms();

	/// Remove the currently selected mesh.
	void removeMesh();

	/// Calls the FileDialog to save a mesh to a file.
	int writeToFile() const;

	/**
	 * checks the mesh quality
	 */
	void checkMeshQuality();

signals:
	void elementSelected(const vtkUnstructuredGridAlgorithm*, bool, int);
	void enableSaveButton(bool);
	void enableRemoveButton(bool);
	void openMeshFile(int);
	void qualityCheckRequested(VtkMeshSource*);
	void removeSelectedMeshComponent();
	void requestCondSetupDialog(const std::string&, const GeoLib::GEOTYPE, const std::size_t, bool on_points);
	void requestMeshRemoval(const QModelIndex&);
	void requestDIRECTSourceTerms(const std::string, const std::vector<GeoLib::Point*>*);
	void saveMeshAction();

/*
    void itemSelectionChanged(const QItemSelection &selected,
        const QItemSelection &deselected);
    //void itemSelectionChangedFromOutside(const QItemSelection &selected,
        const QItemSelection &deselected);
 */
};
#endif // MSHVIEW_H
