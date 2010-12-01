/**
 * \file TreeItem.h
 * KR Initial implementation
 */

#ifndef QTREEITEM_H
#define QTREEITEM_H

#include <QList>
#include <QVariant>

/**
 * \brief Objects nodes for the TreeModel.
 *
 * The TreeItem class provides simple items that contain several pieces of data, 
 * and which can provide information about their parent and child items
 * \sa TreeModel
 */
class TreeItem
{
public:
	TreeItem(const QList<QVariant> &data, TreeItem *parent);
	virtual ~TreeItem();

	void appendChild(TreeItem *child);
	TreeItem* child(int row) const;
	virtual int childCount() const;
	virtual int columnCount() const;
	virtual QVariant data(int column) const;
	virtual bool setData(int column, const QVariant &value);
	int row() const;
	TreeItem* parentItem() const;
	bool removeChildren(int position, int count);

private:
	QList<TreeItem*> _childItems;
	QList<QVariant> _itemData;
	TreeItem* _parentItem;     

};

#endif //QTREEITEM_H
