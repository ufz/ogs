/**
 * \file ColorTableView.h
 * 17/06/2010 KR Initial implementation
 */
#ifndef COLORTABLEVIEW_H
#define COLORTABLEVIEW_H

#include <QTableView>
#include <QItemDelegate>

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
	ColorTableViewDelegate(QWidget *parent = 0) : QItemDelegate(parent) {};

	/// Overwrites the paint-method to set user-defined properties instead of the default properties.
	void paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const;

	QSize sizeHint( const QStyleOptionViewItem &option, const QModelIndex &index ) const;
};

#endif // COLORTABLEVIEW_H

