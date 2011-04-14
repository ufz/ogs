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
class CondObjectListItem;
class vtkPolyDataAlgorithm;

namespace GEOLIB {
	class GeoObject;
	class GeoType;
}

/**
 * \brief The ConditionModel handels FEM conditions such as ICs, BCs and STs on geometric objects.
 */

class ConditionModel : public TreeModel
{
	Q_OBJECT

public:
	ConditionModel(ProjectData &project, QObject* parent = 0);
	~ConditionModel();

    int columnCount(const QModelIndex& parent = QModelIndex()) const;
	vtkPolyDataAlgorithm* vtkSource(const std::string &name, FEMCondition::CondType type);

public slots:
	/// Adds new FEM conditions
	void addConditions(std::vector<FEMCondition*> &conditions);

	/// Removes conditions associated with the given geometry and type.
	void removeFEMConditions(const QString &geometry_name, FEMCondition::CondType type = FEMCondition::UNSPECIFIED);

private:
	/// Adds a new FEM condition to the condition tree model
	void addConditionItem(FEMCondition* conditions);

	/// Removes the FEM condition with the given index.
	//bool removeConditionItem(const QModelIndex &idx);

	CondObjectListItem* createCondParent(TreeItem* parent, FEMCondition::CondType type);
	CondObjectListItem* getCondParent(TreeItem* parent, FEMCondition::CondType type);
	TreeItem* getGEOParent(const QString &geoName, bool create_item = false);
	const GEOLIB::GeoObject* getGEOObject(const std::string &geo_name, GEOLIB::GEOTYPE type, const std::string &obj_name) const;
	size_t getGEOIndex(const std::string &geo_name, GEOLIB::GEOTYPE type, const std::string &obj_name) const;

	ProjectData& _project;

signals:
	void conditionAdded(ConditionModel*, const std::string &name, FEMCondition::CondType);
	void conditionsRemoved(ConditionModel*, const std::string &name, FEMCondition::CondType);
};

#endif // CONDITIONMODEL_H
