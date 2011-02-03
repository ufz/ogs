/**
 * \file PolylinesModel.h
 * 24/9/2009 LB Initial implementation
 * 05/05/2010 KR 2d graphic functionality removed and various layout changes
 */

#ifndef LINESMODEL_H
#define LINESMODEL_H

#include "Model.h"
#include "GeoTreeItem.h"
#include "PolylineVec.h"

#include <QVector>

/**
 * \brief The PolylinesModel is a Qt model which represents geometric Polylines. 
 */
class PolylinesModel : public Model
{
	Q_OBJECT

public:
	PolylinesModel( QString name, const GEOLIB::PolylineVec* polylineVec, QObject* parent = 0 );
	~PolylinesModel();

	int columnCount(const QModelIndex& parent = QModelIndex()) const;

public slots:
	/// Opens a polyline-edit dialog.
	void callEditPlyDialog();
	/// Reloads all items.
	void updateData();

private:
	/// Reconstructs the VtkPolylinesSource object
	void constructVTKObject();
	/// Inserts the polylines data into the TreeModel
	void setData(const GEOLIB::PolylineVec* polylineVec, TreeItem* parent);

	const GEOLIB::PolylineVec* _polylineVec;

private slots:
	/// Calls all necessary functions to connect polyline-segments and update all views and windows.
	void connectPolylineSegments(std::vector<size_t> indexlist, bool, bool);


signals:
	void requestAppendPolylines(const std::vector<GEOLIB::Polyline*> &polylines, std::string &ply_vec_name);
};
#endif // LINESMODEL_H
