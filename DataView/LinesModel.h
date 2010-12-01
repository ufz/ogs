/**
 * \file PolylinesModel.h
 * 24/9/2009 LB Initial implementation
 * 05/05/2010 KR 2d graphic functionality removed and various layout changes
 */

#ifndef LINESMODEL_H
#define LINESMODEL_H

#include "Model.h"
#include "TreeItem.h"
#include "PolylineVec.h"

#include <QVector>

/**
 * The PolylinesModel is a Qt model which represents Polylines. 
 */
class PolylinesModel : public Model
{
	Q_OBJECT

public:
	PolylinesModel( QString name, const GEOLIB::PolylineVec* polylineVec, QObject* parent = 0 );
	~PolylinesModel();

	int columnCount(const QModelIndex& parent = QModelIndex()) const;

public slots:
	/// Reloads all items.
	void updateData();

private:
	void setData(const GEOLIB::PolylineVec* polylineVec, TreeItem* parent);

	const GEOLIB::PolylineVec* _polylineVec;

};
#endif // LINESMODEL_H
