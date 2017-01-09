/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the TreeModel class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QAbstractItemModel>

class QVariant;
class QModelIndex;
class TreeItem;

/**
 * \brief A hierarchical model for a tree implemented as a double-linked list
 *
 * A hierarchical model for the pre-defined QTreeView included in QT. The tree as implemented
 * as a double-linked list.
 */
class TreeModel : public QAbstractItemModel
{
    Q_OBJECT

public:
    TreeModel(QObject* parent = 0);
    virtual ~TreeModel();

    QVariant data(const QModelIndex &index, int role) const;
    bool setData(const QModelIndex &index, const QVariant &value, int role /* = Qt::EditRole */);
    Qt::ItemFlags flags(const QModelIndex &index) const;
    TreeItem* getItem(const QModelIndex &index) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role =
                                Qt::DisplayRole) const;
    QModelIndex index(int row, int column, const QModelIndex &parent = QModelIndex()) const;
    QModelIndex parent(const QModelIndex &index) const;
    bool removeRows(int row, int count, const QModelIndex & parent);
    int rowCount(const QModelIndex &parent = QModelIndex()) const;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;

    TreeItem* rootItem() const;

public slots:
    void updateData();

protected:
    TreeItem* _rootItem;

private:
    void setupModelData(const QStringList &lines, TreeItem* parent);
};
