/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file DataView.h
 *
 * Created on 2009-09-24 by Lars Bilke
 */
#ifndef DATAVIEW_H
#define DATAVIEW_H

#include "Point.h"
#include "GeoType.h"
#include <QTreeView>

class MshModel;
class VtkMeshSource;

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
class DataView : public QTreeView
{
	Q_OBJECT

public:
	DataView(QWidget* parent = 0);
	~DataView();

public slots:
	void updateView();

protected slots:
	/// Is called when the selection of this view changes. Emits a the signal
	/// itemSelectionChanged()
	//void selectionChanged(const QItemSelection &selected, const QItemSelection &deselected);

	/// Selects items without sending signals.
	//void selectionChangedFromOutside(const QItemSelection &selected,
	//	const QItemSelection &deselected);

	/// Clears the selection
	//void clearSelection();

private:
	void contextMenuEvent( QContextMenuEvent* event );

private slots:
	/// Open a dialog for editing meshes.
	void openMshEditDialog();

	/// Adds a new mesh.
	void addMeshAction();

	void addDIRECTSourceTerms();

	void loadDIRECTSourceTerms();

	/// Remove the currently selected mesh.
	void removeMesh();

	/// Remove all currently loaded meshes.
	void removeAllMeshes();

	/// Calls the FileDialog to save a mesh to a file.
	int writeMeshToFile() const;

	/**
	 * checks the mesh quality
	 */
	void checkMeshQuality();

signals:
	void qualityCheckRequested(VtkMeshSource*);
	void requestCondSetupDialog(const std::string&, const GeoLib::GEOTYPE, const size_t, bool on_points);
	void requestMeshRemoval(const QModelIndex&);
	void requestDIRECTSourceTerms(const std::string, const std::vector<GeoLib::Point*>*);
	void saveMeshAction();

/*
    void itemSelectionChanged(const QItemSelection &selected,
        const QItemSelection &deselected);
    void itemSelectionChangedFromOutside(const QItemSelection &selected,
        const QItemSelection &deselected);
 */
};
#endif // DATAVIEW_H
