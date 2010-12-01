/**
 * \file ConditionModel.h
 * 18/10/2010 KR Initial implementation
 */

#ifndef CONDITIONMODEL_H
#define CONDITIONMODEL_H

// ** INCLUDES **
#include "TreeModel.h"
#include "ProjectData.h"
#include <QVector>

class CSourceTerm;

/**
 * The ConditionModel handels conditions such as ICs, BCs and STs on geometric objects 
 */

class ConditionModel : public TreeModel
{
	Q_OBJECT

public:
	ConditionModel(ProjectData &project, QObject* parent = 0);
	~ConditionModel();

    int columnCount(const QModelIndex& parent = QModelIndex()) const;
	//QVariant data(const QModelIndex& index, int role) const;
	//QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;

	//bool setData(const QModelIndex& index, const QVariant& value, int role = Qt::EditRole);

public slots:
	/// Adds new source terms
	void addSourceTerms(std::vector<CSourceTerm*> *sourceterms, const QString &name);
	/// Returns the mesh with the given index.
	const std::vector<CSourceTerm*> *getSourceTerms(const QModelIndex &idx) const;
	/// Reloads all items.
	//void updateData();

private:
	//void setData(std::vector<GEOLIB::GeoObject*> *points, TreeItem* parent);

	ProjectData& _project;
	TreeItem* _bcParent;
	TreeItem* _icParent;
	TreeItem* _stParent;

};

#endif // CONDITIONMODEL_H
