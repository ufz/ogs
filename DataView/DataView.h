/**
 * \file DataView.h
 * 24/9/2009 LB Initial implementation
 */
#ifndef DATAVIEW_H
#define DATAVIEW_H

#include <QTreeView>

class MshModel;
class VtkMeshSource;

namespace Mesh_Group {
	class CFEMesh;
}

/**
 *	The DataView is table view which acts as a base class for displaying
 *  several OSG data formats.
 */
class DataView : public QTreeView
{
	Q_OBJECT

public:
	DataView(QWidget* parent = 0);

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
	void requestMeshRemoval(const QModelIndex&);
	void saveMeshAction();
/*
	void itemSelectionChanged(const QItemSelection &selected,
		const QItemSelection &deselected);
	void itemSelectionChangedFromOutside(const QItemSelection &selected,
		const QItemSelection &deselected);
*/
};
#endif // DATAVIEW_H
