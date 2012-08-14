/**
 * \file ColorTableModel.h
 * 17/06/2010 KR Initial implementation
 */

#ifndef COLORTABLEMODEL_H
#define COLORTABLEMODEL_H

#include "Color.h"
#include <QAbstractTableModel>
#include <QColor>

/**
 * The PolylinesModel is a Qt model which represents Polylines.
 */
class ColorTableModel : public QAbstractTableModel
{
	Q_OBJECT

public:
	ColorTableModel( const std::map<std::string, GeoLib::Color*> &colorLookupTable,
	                 QObject* parent = 0 );
	~ColorTableModel();

	int columnCount(const QModelIndex& parent = QModelIndex()) const;

	QVariant data( const QModelIndex& index, int role ) const;

	int rowCount(const QModelIndex& parent = QModelIndex()) const
	{
		Q_UNUSED (parent);
		return _listOfPairs.size();
	}

	QVariant headerData( int section, Qt::Orientation orientation,
	                     int role /*= Qt::DisplayRole*/ ) const;

private:
	bool buildTable( const std::map<std::string, GeoLib::Color*> &colorLookupTable );

	QList< QPair<QString, QColor> > _listOfPairs;
};
#endif // COLORTABLEMODEL_H
