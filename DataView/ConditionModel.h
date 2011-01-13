/**
 * \file ConditionModel.h
 * 18/10/2010 KR Initial implementation
 */

#ifndef CONDITIONMODEL_H
#define CONDITIONMODEL_H

// ** INCLUDES **
#include "TreeModel.h"
#include "ProjectData.h"

class FEMCondition;

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
	/// Adds new FEM condition
	void addCondition(FEMCondition* conditions);

	/// Adds new FEM conditions
	void addConditions(std::vector<FEMCondition*> &conditions);

	/// Returns the FEM condition set on a GeoObject with the given name.
	//const FEMCondition* getCondition(const string &name) const;

	/// Removes the FEM condition set on a GeoObject with the given name.
	//bool removeCondition(const string &name);

	/// Removes the FEM condition with the given index.
	bool removeCondition(const QModelIndex &idx);

	/// Reloads all items.
	//void updateData();

private:
	//void setData(std::vector<GEOLIB::GeoObject*> *points, TreeItem* parent);

	ProjectData& _project;
	TreeItem* _bcParent;
	TreeItem* _icParent;
	TreeItem* _stParent;

signals:
	void condAdded(ConditionModel*, const QModelIndex&);
	void condRemoved(ConditionModel*, const QModelIndex&);
};

#endif // CONDITIONMODEL_H
