/**
 * \file SurfaceModel.h
 * 23/04/2010 KR Initial implementation
 *
 */

#ifndef SURFACEMODEL_H
#define SURFACEMODEL_H

#include "Model.h"
#include "TreeItem.h"
#include "SurfaceVec.h"

#include <QVector>

/**
 * The SurfaceModel is a Qt model which represents Polylines. 
 */
class SurfaceModel : public Model
{
	Q_OBJECT

public:
	SurfaceModel( QString name, const GEOLIB::SurfaceVec* surfaceVec, QObject* parent = 0 );
	~SurfaceModel();

	int columnCount(const QModelIndex& parent = QModelIndex()) const;

public slots:
	/// Reloads all items.
	void updateData();

private:
	void setData(const GEOLIB::SurfaceVec* surfaceVec, TreeItem* parent);

	const GEOLIB::SurfaceVec* _surfaceVec;

};
#endif // SURFACEMODEL_H
