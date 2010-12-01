/**
 * \file ColorTableModel.cpp
 * 24/9/2009 LB Initial implementation
 * 05/05/2010 KR 2d graphic functionality removed and various layout changes
 *
 * Implementation of PolylinesModel
 */

#include "ColorTableModel.h"


ColorTableModel::ColorTableModel( const std::map<std::string, GEOLIB::Color*> &colorLookupTable, QObject* parent /*= 0*/ )
{
	Q_UNUSED(parent)

	this->buildTable(colorLookupTable);
}

ColorTableModel::~ColorTableModel()
{
}

int ColorTableModel::columnCount( const QModelIndex& parent /*= QModelIndex()*/ ) const
{
	Q_UNUSED(parent)

	return 2;
}

QVariant ColorTableModel::headerData( int section, Qt::Orientation orientation, int role /*= Qt::DisplayRole*/ ) const
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

    if (index.row() >= _listOfPairs.size() || index.row()<0)
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

bool ColorTableModel::buildTable(const std::map<std::string, GEOLIB::Color*> &colorLookupTable)
{
     int count = 0;
     beginInsertRows(QModelIndex(), 0, colorLookupTable.size()-1);

	 for (std::map<std::string, GEOLIB::Color*>::const_iterator it=colorLookupTable.begin(); it !=colorLookupTable.end(); ++it)
	 {
		 QColor color((*(it->second))[0], (*(it->second))[1], (*(it->second))[2]);
		 QString name(QString::fromStdString(it->first));

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


