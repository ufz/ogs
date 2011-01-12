/**
 * \file ConditionView.h
 * 2010/12/13 KR Initial implementation
 */

#ifndef CONDITIONVIEW_H
#define CONDITIONVIEW_H

#include <QTreeView>
#include <QContextMenuEvent>

/**
 * \brief A view for the StationTreeModel with a number of properties adequate for this kind of data
 * \sa StationTreeModel, ModelTreeItem
 */
class ConditionView : public QTreeView
{
	Q_OBJECT

public:
	/// Constructor
	ConditionView(QWidget* parent = 0);

	/// Update the view to visualise changes made to the underlying data
	void updateView();

protected slots:
	/// Instructions if the selection of items in the view has changed.
	void selectionChanged(const QItemSelection &selected, const QItemSelection &deselected);

private:
	/// Actions to be taken after a right mouse click is performed in the station view.
	void contextMenuEvent( QContextMenuEvent* e );

private slots:
	void on_Clicked(QModelIndex idx);
	void removeCondition();

signals:
	void itemSelectionChanged(const QItemSelection & selected, const QItemSelection & deselected);
	void conditionRemoved(std::string name);
};

#endif //CONDITIONVIEW_H
