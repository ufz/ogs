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
class ProcessItem;
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
	/// Returns the vtk source object for the specified subtree of a process with the given name.
	vtkPolyDataAlgorithm* vtkSource(const FiniteElement::ProcessType pcs_type, const FEMCondition::CondType cond_type);

public slots:
	/// Adds a vector of FEM Conditions to the model. Objects in the vector can consist of BCs, ICs or STs in any combination and sequence.
	void addConditions(std::vector<FEMCondition*> &conditions);

	/// Adds a single FEM Conditions to the model
	void addCondition(FEMCondition* condition);

	/// Adds a process to the model
	ProcessItem* addProcess(ProcessInfo* pcs);

	/// Removes FEMConditions from the the model. Conditions can be specified by process type, geometry name or condition type or a combination of the three.
	void removeFEMConditions(const FiniteElement::ProcessType pcs_type, const std::string &geometry_name, const FEMCondition::CondType cond_type);

	/// Removes a process from the model
	void removeProcess(const FiniteElement::ProcessType type);

	/// Removes all processes from the model
	void removeAllProcesses();

private:
	/// Adds a new FEM condition to the condition tree model.
	void addConditionItem(FEMCondition* condition);

	/// Removes the FEM condition with the given index.
	//bool removeConditionItem(const QModelIndex &idx);

	/// Creates the TreeItem for one of the condition subtrees.
	CondObjectListItem* createCondParent(ProcessItem* parent, const FEMCondition::CondType type, const std::string &geometry_name);

	/// Returns the subtree-item for a given type of condtion.
	CondObjectListItem* getCondParent(TreeItem* parent, const FEMCondition::CondType type) ;

	/// Returns the subtree item for a process with the given name. If create_item is true, this item will be created if it doesn't exist yet.
	ProcessItem* getProcessParent(const FiniteElement::ProcessType type) const;

	/// Returns the index of a geometric item of the given name and type for the associated geometry.
	int getGEOIndex(const std::string &geo_name,
	                GEOLIB::GEOTYPE type,
	                const std::string &obj_name) const;

	ProjectData& _project;

signals:
	void conditionAdded(ProcessModel*, const FiniteElement::ProcessType, const FEMCondition::CondType);
	void conditionsRemoved(ProcessModel*, const FiniteElement::ProcessType, const FEMCondition::CondType);
};

#endif // PROCESSMODEL_H
