/**
 * \file ProcessModel.h
 * 18/10/2010 KR Initial implementation
 */

#ifndef PROCESSMODEL_H
#define PROCESSMODEL_H

// ** INCLUDES **
#include "ProjectData.h"
#include "TreeModel.h"

class FEMCondition;
class CondObjectListItem;
class vtkPolyDataAlgorithm;

namespace GEOLIB
{
class GeoObject;
}

/**
 * \brief A model implementing a tree structure for process-relevant information such as
 * process types, FEM-Conditions (BCs, ICs, STs), etc. as a double-linked list.
 * \sa TreeModel, ProcessView, TreeItem, CondObjectListItem
 */
class ProcessModel : public TreeModel
{
	Q_OBJECT

public:
	ProcessModel(ProjectData &project, QObject* parent = 0);
	~ProcessModel();

	int columnCount(const QModelIndex& parent = QModelIndex()) const;
	/// Returns the vtk source object for the specified subtree of a geometry with the given name.
	vtkPolyDataAlgorithm* vtkSource(const std::string &name, FEMCondition::CondType type);

public slots:
	/// Adds a vector of FEM Conditions to the model. Objects in the vector can consist of BCs, ICs or STs in any combination and sequence.
	void addConditions(std::vector<FEMCondition*> &conditions);

	/// Adds a single FEM Conditions to the model
	void addCondition(FEMCondition* condition);

	/// Removes a subtree (BCs, ICs, STs) from the the model. If all conditions for a given geometry are removed, this tree is completely removed.
	void removeFEMConditions(const QString &geometry_name,
	                         FEMCondition::CondType type = FEMCondition::UNSPECIFIED);

private:
	/// Adds a new FEM condition to the condition tree model.
	void addConditionItem(FEMCondition* condition);

	/// Removes the FEM condition with the given index.
	//bool removeConditionItem(const QModelIndex &idx);

	/// Creates the TreeItem for one of the condition subtrees.
	CondObjectListItem* createCondParent(TreeItem* parent, FEMCondition::CondType type);

	/// Returns the subtree-item for a given type of condtion.
	CondObjectListItem* getCondParent(TreeItem* parent, FEMCondition::CondType type);

	/// Returns the subtree item for a geometry with the given name. If create_item is true, this item will be created if it doesn't exist yet.
	TreeItem* getGEOParent(const QString &geoName, bool create_item = false);

	/// Returns the index of a geometric item of the given name and type for the associated geometry.
	int getGEOIndex(const std::string &geo_name,
	                GEOLIB::GEOTYPE type,
	                const std::string &obj_name) const;

	ProjectData& _project;

signals:
	void conditionAdded(ProcessModel*, const std::string &name, FEMCondition::CondType);
	void conditionsRemoved(ProcessModel*, const std::string &name, FEMCondition::CondType);
};

#endif // PROCESSMODEL_H
