/**
 * \file ConditionView.h
 * 2010/12/13 KR Initial implementation
 */

#ifndef CONDITIONVIEW_H
#define CONDITIONVIEW_H

#include <QContextMenuEvent>
#include <QTreeView>

#include "FEMCondition.h"

class ConditionModel;

/**
 * \brief A view for FEM-Conditions (Initial- & Boundary Conditions / Source Terms) with a number of additional
 * information such as Process Type, Distribution, etc.
 * \sa ConditionModel, CondItem
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
	//void removeAllConditions();

signals:
	void itemSelectionChanged(const QItemSelection & selected,
	                          const QItemSelection & deselected);
	void conditionsRemoved(QString, FEMCondition::CondType);
};

#endif //CONDITIONVIEW_H
