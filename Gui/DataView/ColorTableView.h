/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file ColorTableView.h
 *
 * Created on 2010-06-17 by Karsten Rink
 */
#ifndef COLORTABLEVIEW_H
#define COLORTABLEVIEW_H

#include <QItemDelegate>
#include <QTableView>

/**
 *	A QTableView to display colour lookup tables.
 */
class ColorTableView : public QTableView
{
	Q_OBJECT

public:
	/// Constructor
	ColorTableView(QWidget* parent = 0);
};

/**
 *	A delegate class to manage properties of ColorTableView.
 */
class ColorTableViewDelegate : public QItemDelegate
{
	Q_OBJECT

public:
	/// Constructor
	ColorTableViewDelegate(QWidget* parent = 0) : QItemDelegate(parent) {}

	/// Overwrites the paint-method to set user-defined properties instead of the default properties.
	void paint(QPainter* painter, const QStyleOptionViewItem &option,
	           const QModelIndex &index) const;

	QSize sizeHint( const QStyleOptionViewItem &option, const QModelIndex &index ) const;
};

#endif // COLORTABLEVIEW_H

