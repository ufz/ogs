/**
 * \file
 * \author Lars Bilke
 * \date   2009-09-24
 * \brief  Implementation of the ColorTableModel class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ColorTableModel.h"

ColorTableModel::ColorTableModel( const std::map<std::string, DataHolderLib::Color*> &colorLookupTable,
                                  QObject* parent /*= 0*/ )
{
    Q_UNUSED(parent)

    this->buildTable(colorLookupTable);
}

ColorTableModel::~ColorTableModel() = default;

int ColorTableModel::columnCount( const QModelIndex& parent /*= QModelIndex()*/ ) const
{
    Q_UNUSED(parent)

    return 2;
}

QVariant ColorTableModel::headerData( int section, Qt::Orientation orientation,
                                      int role /*= Qt::DisplayRole*/ ) const
{
    if (role != Qt::DisplayRole)
        return QVariant();

    if (orientation == Qt::Horizontal)
    {
        switch (section)
        {
        case 0: return "Name";
        case 1: return "Colour";
        default: return QVariant();
        }
    }
    else
        return QString("Row %1").arg(section);
}

QVariant ColorTableModel::data( const QModelIndex& index, int role ) const
{
    if (!index.isValid())
        return QVariant();

    if (index.row() >= _listOfPairs.size() || index.row() < 0)
        return QVariant();

    if (role == Qt::DisplayRole)
    {
        QPair<QString, QColor> pair = _listOfPairs.at(index.row());

        switch (index.column())
        {
        case 0:
            return pair.first;
        case 1:
            return pair.second;
        default:
            return QVariant();
        }
    }
    return QVariant();
}

bool ColorTableModel::buildTable(const std::map<std::string, DataHolderLib::Color*> &colorLookupTable)
{
    int count = 0;
    beginInsertRows(QModelIndex(), 0, colorLookupTable.size() - 1);

    for (const auto& row : colorLookupTable)
    {
        QColor color((*row.second)[0], (*row.second)[1], (*row.second)[2]);
        QString name(QString::fromStdString(row.first));

        /* Saudi Arabia strat names *
           if (it->first.compare("1")==0) name="Buweib";
           if (it->first.compare("2")==0) name="Wasia";
           if (it->first.compare("3")==0) name="Aruma";
           if (it->first.compare("4")==0) name="Umm Er Radhuma";
           if (it->first.compare("5")==0) name="Rus";
           if (it->first.compare("6")==0) name="Dammam";
           if (it->first.compare("7")==0) name="Neogene";
         */

        QPair<QString, QColor> pair(name, color);
        _listOfPairs.insert(count++, pair);
    }

    endInsertRows();
    return true;
}

