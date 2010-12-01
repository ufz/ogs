/**
 * \file DataView.h
 * 24/9/2009 LB Initial implementation
 */
#ifndef DATAVIEW_H
#define DATAVIEW_H

#include <QTreeView>

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

	/// Returns the selected indexes. Overwritten from QTableView to make it public.
	//QModelIndexList selectedIndexes() const;

protected slots:
	/// Is called when the selection of this view changes. Emits a the signal
	/// itemSelectionChanged()
	//void selectionChanged(const QItemSelection &selected,
	//	const QItemSelection &deselected);

	/// Selects items without sending signals.
	//void selectionChangedFromOutside(const QItemSelection &selected,
	//	const QItemSelection &deselected);

	/// Clears the selection
	//void clearSelection();

private:

signals:
	void itemSelectionChanged(const QItemSelection &selected,
		const QItemSelection &deselected);
	void itemSelectionChangedFromOutside(const QItemSelection &selected,
		const QItemSelection &deselected);
	
};
#endif // DATAVIEW_H
