/**
 * \file CheckboxDelegate.h
 * 19/08/2010 LB Initial implementation
 */

#ifndef CHECKBOXDELEGATE_H
#define CHECKBOXDELEGATE_H

#include <QItemDelegate>

class QWidget;
class QRect;

/**
 * \brief CheckboxDelegate modifies a model view to display boolean values as checkboxes.
 *
 * Important: the column on which this delegate is set (QAbstractItemView::setItemDelegateForColumn())
 * must not have the flags Qt::ItemIsEditable or Qt::ItemIsUserCheckable set in the model.
 **/
class CheckboxDelegate : public QItemDelegate
{
	Q_OBJECT

public:
	/// \brief Constructor
	CheckboxDelegate (QObject* parent = 0);

	/// \brief Paints a checkbox. This overrides the default painting of a combo box.
	void paint(QPainter* painter, const QStyleOptionViewItem& option,
	           const QModelIndex& index) const;

	/// \brief Handles the click events and sets the model data.
	bool editorEvent(QEvent* event, QAbstractItemModel* model,
	                 const QStyleOptionViewItem &option, const QModelIndex &index);

	QSize sizeHint (const QStyleOptionViewItem & option, const QModelIndex & index) const;

private:
	QRect checkboxRect(const QStyleOptionViewItem& viewItemStyleOptions) const;
};

#endif // CHECKBOXDELEGATE_H
