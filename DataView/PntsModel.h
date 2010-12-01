/**
 * \file PntsModel.h
 * 24/9/2009 LB Initial implementation
 * 05/05/2010 KR 2d graphic functionality removed and various layout changes
 */

#ifndef PNTSMODEL_H
#define PNTSMODEL_H

// ** INCLUDES **
#include "PointVec.h"
#include "Model.h"
#include <QVector>

/**
 * The PntsModel is a Qt model which represents Point objects.
 */
class PntsModel : public Model
{
	Q_OBJECT

public:
	PntsModel(QString name, const GEOLIB::PointVec* pntVec, QObject* parent = 0);
	~PntsModel();

    int columnCount(const QModelIndex& parent = QModelIndex()) const;
	QVariant data(const QModelIndex& index, int role) const;
	QVariant headerData(int section, Qt::Orientation orientation,
		int role = Qt::DisplayRole) const;

	bool setData(const QModelIndex& index, const QVariant& value,
		int role = Qt::EditRole);

public slots:
	/// Reloads all items.
	void updateData();

private:
	void setData(const GEOLIB::PointVec* pntVec, TreeItem* parent);

	const GEOLIB::PointVec* _pntVec;

};

#endif // PNTSMODEL_H
