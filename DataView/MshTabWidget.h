/**
 * \file MshTabWidget.h
 * 3/11/2009 LB Initial implementation
 * 18/05/2010 KR Re-Implementation
 *
 */


#ifndef MSHTABWIDGET_H
#define MSHTABWIDGET_H

// ** INCLUDES **
#include "ui_MshTabWidgetBase.h"

class MshModel;
class VtkMeshSource;

namespace Mesh_Group {
	class CFEMesh;
}

/**
 * \brief Widget for data views of meshes.
 */
class MshTabWidget : public QWidget, public Ui_MshTabWidgetBase
{
	Q_OBJECT

public:
	MshTabWidget(QWidget* parent = 0);

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

};

#endif // MSHTABWIDGET_H
