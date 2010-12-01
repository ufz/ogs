/**
 * \file QOgsDataView.h
 * 24/9/2009 LB Initial implementation
 */
#ifndef QOGSDATAVIEW_H
#define QOGSDATAVIEW_H

#include <QTableView>

/**
 *	The QOgsDataView is table view which displays CGLPolyline data.
 */
class QOgsDataView : public QTableView
{
	Q_OBJECT

public:
	QOgsDataView(QWidget* parent = 0);

	/// Returns the selected indexes. Overwritten from QTableView to make it public.
	QModelIndexList selectedIndexes() const;

protected slots:
	void selectionChanged(const QItemSelection &selected, const QItemSelection &deselected);
	void selectionChangedFromOutside(const QItemSelection &selected, const QItemSelection &deselected);

private:

signals:
	void itemSelectionChanged(const QItemSelection &selected, const QItemSelection &deselected);
	void itemSelectionChangedFromOutside(const QItemSelection &selected, const QItemSelection &deselected);
	
};
#endif // QOGSDATAVIEW_H
