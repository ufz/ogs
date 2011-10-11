/**
 * \file ConditionModel.h
 * 18/10/2010 KR Initial implementation
 */

#ifndef CONDITIONMODEL_H
#define CONDITIONMODEL_H

// ** INCLUDES **
#include "ProjectData.h"
#include "TreeModel.h"

class FEMCondition;
class CondObjectListItem;
class vtkPolyDataAlgorithm;

namespace GEOLIB
{
class GeoObject;
class GeoType;
}

/**
 * \brief A model for the ConditionView implementing a tree of FEM-Conditions (BCs, ICs, STs) as a double-linked list.
 * \sa TreeModel, ConditionView, TreeItem, CondObjectListItem
 */
class ConditionModel : public TreeModel
{
	Q_OBJECT

public:
	ConditionModel(ProjectData &project, QObject* parent = 0);
	~ConditionModel();

	int columnCount(const QModelIndex& parent = QModelIndex()) const;
	/// Returns the vtk source object for the specified subtree of a geometry with the given name.
	vtkPolyDataAlgorithm* vtkSource(const std::string &name, FEMCondition::CondType type);

public slots:
	/// Adds a vector of FEM Conditions to the model. Objects in the vector can consist of BCs, ICs or STs in any combination and sequence.
	void addConditions(std::vector<FEMCondition*> &conditions);

	/// Removes a subtree (BCs, ICs, STs) from the the model. If all conditions for a given geometry are removed, this tree is completely removed.
	void removeFEMConditions(const QString &geometry_name,
	                         FEMCondition::CondType type = FEMCondition::UNSPECIFIED);

private:
	/// Adds a new FEM condition to the condition tree model.
	void addConditionItem(FEMCondition* conditions);

	/// Removes the FEM condition with the given index.
	//bool removeConditionItem(const QModelIndex &idx);

	/// Creates the TreeItem for one of the condition subtrees.
	CondObjectListItem* createCondParent(TreeItem* parent, FEMCondition::CondType type);

	/// Returns the subtree-item for a given type of condtion.
	CondObjectListItem* getCondParent(TreeItem* parent, FEMCondition::CondType type);

	/// Returns the subtree item for a geometry with the given name. If create_item is true, this item will be created if it doesn't exist yet.
	TreeItem* getGEOParent(const QString &geoName, bool create_item = false);
	/// Returns the geo object for a geometric item of the given name and type for the associated geometry.
	const GEOLIB::GeoObject* getGEOObject(const std::string &geo_name,
	                                      GEOLIB::GEOTYPE type,
	                                      const std::string &obj_name) const;
	/// Returns the index of a geometric item of the given name and type for the associated geometry.
	int getGEOIndex(const std::string &geo_name,
	                GEOLIB::GEOTYPE type,
	                const std::string &obj_name) const;

	ProjectData& _project;

signals:
	void conditionAdded(ConditionModel*, const std::string &name, FEMCondition::CondType);
	void conditionsRemoved(ConditionModel*, const std::string &name, FEMCondition::CondType);
};

#endif // CONDITIONMODEL_H
