/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-07
 * \brief  Implementation of the FEMConditionSetupDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FEMConditionSetupDialog.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "OGSError.h"
#include "FEMEnums.h"
#include "ProjectData.h"
#include "StrictDoubleValidator.h"
#include "LinearEditDialog.h"
#include "CondFromRasterDialog.h"

#include "BoundaryCondition.h"
#include "InitialCondition.h"
#include "SourceTerm.h"


FEMConditionSetupDialog::FEMConditionSetupDialog(const std::string &associated_geometry,
												 const GeoLib::GEOTYPE type,
												 const std::string &geo_name,
												 const GeoLib::GeoObject* const geo_object,
												 bool  on_points,
												 QDialog* parent)
: QDialog(parent), _cond(associated_geometry, FEMCondition::UNSPECIFIED), _set_on_points(on_points),
  _combobox(nullptr), directButton(nullptr), _mesh(nullptr), _first_value_validator(nullptr)
{
	_cond.setGeoType(type);
	_cond.setGeoName(geo_name);
	_cond.setGeoObj(geo_object);

	setupUi(this);
	setupDialog();
}

FEMConditionSetupDialog::FEMConditionSetupDialog(const FEMCondition &cond, QDialog* parent)
	: QDialog(parent), _cond(cond), _set_on_points(false), _combobox(nullptr), directButton(nullptr),
	_mesh(nullptr), _first_value_validator(nullptr)
{
	setupUi(this);
	setupDialog();
	setValuesFromCond();
}

FEMConditionSetupDialog::FEMConditionSetupDialog(const std::string &name, const MeshLib::Mesh* mesh, QDialog* parent)
: QDialog(parent), _cond(name, FEMCondition::UNSPECIFIED), _set_on_points(false),  _combobox(nullptr), directButton(nullptr),
  _mesh(mesh), _first_value_validator(nullptr)
{
	_cond.setGeoType(GeoLib::GEOTYPE::INVALID);
	_cond.setGeoName(name);
	_cond.setGeoObj(nullptr);

	setupUi(this);
	setupDialog();
}

FEMConditionSetupDialog::~FEMConditionSetupDialog()
{
	delete _combobox;
	delete directButton;
	delete _first_value_validator;
}

void FEMConditionSetupDialog::setupDialog()
{
	if (_cond.getGeomType() != GeoLib::GEOTYPE::INVALID)
	{
		this->disTypeBox->addItem("Constant (Dirichlet)");
		if (_cond.getGeomType() == GeoLib::GEOTYPE::POLYLINE)
			this->disTypeBox->addItem("Linear (Dirichlet)");

		if (this->_set_on_points)
		{
			_combobox = new QComboBox;
			_combobox->addItem("Elevation");
			static_cast<QGridLayout*>(this->layout())->addWidget(_combobox,5,1) ;
		}
		else
		{
			_first_value_validator = new StrictDoubleValidator(-1e+10, 1e+10, 5);
			this->firstValueEdit->setText("0");
			this->firstValueEdit->setValidator (_first_value_validator);
		}
	}
	else	// direct on mesh
	{
		this->disTypeBox->addItem("Direct");
		this->setValueInputWidget(true);
		this->condTypeBox->setCurrentIndex(2);
	}

	const std::list<std::string> process_names = FiniteElement::getAllProcessNames();
	for (std::list<std::string>::const_iterator it = process_names.begin(); it != process_names.end(); ++it)
		this->processTypeBox->addItem(QString::fromStdString(*it));

	const std::list<std::string> pv_names = FiniteElement::getAllPrimaryVariableNames();
	for (std::list<std::string>::const_iterator it = pv_names.begin(); it != pv_names.end(); ++it)
		this->pvTypeBox->addItem(QString::fromStdString(*it));
/*
	const std::list<std::string> dis_names = FiniteElement::getAllDistributionNames();
	for (std::list<std::string>::const_iterator it = dis_names.begin(); it != dis_names.end(); ++it)
		this->disTypeBox->addItem(QString::fromStdString(*it));
*/
}

void FEMConditionSetupDialog::setValuesFromCond()
{
	QString pcs_type = QString::fromStdString(FiniteElement::convertProcessTypeToString(_cond.getProcessType()));
	this->processTypeBox->setCurrentIndex(this->processTypeBox->findText(pcs_type));

	QString pv_type = QString::fromStdString(FiniteElement::convertPrimaryVariableToString(_cond.getProcessPrimaryVariable()));
	this->pvTypeBox->setCurrentIndex(this->pvTypeBox->findText(pv_type));

	if (_cond.getCondType() == FEMCondition::INITIAL_CONDITION)
		this->condTypeBox->setCurrentIndex(1);
	else if (_cond.getCondType() == FEMCondition::SOURCE_TERM)
	{
		this->condTypeBox->setCurrentIndex(2);
		on_condTypeBox_currentIndexChanged(2);
	}

	if (_cond.getGeomType() != GeoLib::GEOTYPE::INVALID)
	{
		if (_cond.getProcessDistributionType() == FiniteElement::CONSTANT || _cond.getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN)
		{
			this->disTypeBox->setCurrentIndex(0);
			this->firstValueEdit->setText(QString::number(_cond.getDisValues()[0]));
		}
		else
		{
			this->disTypeBox->setCurrentIndex(1);
			directButton = new QPushButton(QString::number(static_cast<int>(_cond.getDisValues().size())) + " values");
		}
	}
	else	// direct on mesh
	{
		this->directButton->setText(QString::number(static_cast<int>(_cond.getDisValues().size())) + " values");
	}
}


void FEMConditionSetupDialog::accept()
{
	_cond.setProcessType(static_cast<FiniteElement::ProcessType>(this->processTypeBox->currentIndex() + 1));
	_cond.setProcessPrimaryVariable(static_cast<FiniteElement::PrimaryVariable>(this->pvTypeBox->currentIndex() + 1));

	if (_cond.getGeomType() != GeoLib::GEOTYPE::INVALID)
	{
		if (condTypeBox->currentIndex()>1) // ST
		{
			if (this->disTypeBox->currentIndex()>0)
				_cond.setProcessDistributionType(FiniteElement::LINEAR_NEUMANN);
			else
			{
				_cond.setProcessDistributionType(FiniteElement::CONSTANT_NEUMANN);
				_cond.setConstantDisValue(this->firstValueEdit->text().toDouble());
			}
		}
		else // BC or IC
		{
			if (this->disTypeBox->currentIndex()>0)
				_cond.setProcessDistributionType(FiniteElement::LINEAR);
			else
			{
				_cond.setProcessDistributionType(FiniteElement::CONSTANT);
				_cond.setConstantDisValue(this->firstValueEdit->text().toDouble());
			}
		}
	}
	else	// direct on mesh
	{
		if (this->condTypeBox->currentIndex()==1) // IC
		{
			_cond.setProcessDistributionType(FiniteElement::NODESCONSTANT);
			_cond.setConstantDisValue(this->firstValueEdit->text().toDouble());
		}
		else // BC or ST
			_cond.setProcessDistributionType(FiniteElement::DIRECT);
	}
	if (_cond.getDisValues().size()==0)
	{
		OGSError::box("No distribution values specified!");
		return;
	}

	if (!_set_on_points)
	{
		std::vector<FEMCondition*> conditions;
		conditions.push_back(this->typeCast(_cond));
		emit createFEMCondition(conditions);
	}
	else
		this->copyCondOnPoints();

	this->done(QDialog::Accepted);
}

void FEMConditionSetupDialog::reject()
{
	this->done(QDialog::Rejected);
}


void FEMConditionSetupDialog::on_condTypeBox_currentIndexChanged(int index)
{
	//if (index==1)
	//	this->geoNameBox->addItem("Domain");
	// remove "Domain" if IC is unselected
	if (_cond.getGeomType() != GeoLib::GEOTYPE::INVALID)
	{
		if (index>1) // source terms selected
		{
			this->clearDisTypeBox();
			this->disTypeBox->addItem("Constant (Neumann)");
			if (_cond.getGeomType() == GeoLib::GEOTYPE::POLYLINE)
				this->disTypeBox->addItem("Linear (Neumann)");
		}
		else
		{
			this->clearDisTypeBox();
			this->disTypeBox->addItem("Constant (Dirichlet)");
			if (_cond.getGeomType() == GeoLib::GEOTYPE::POLYLINE)
				this->disTypeBox->addItem("Linear (Dirichlet)");
		}
	}
	else {
		if (index==1) // initial condition selected
		{
			this->clearDisTypeBox();
			this->disTypeBox->addItem("Domain");
			this->setValueInputWidget(false);
		}
		else
		{
			this->clearDisTypeBox();
			this->disTypeBox->addItem("Direct");
			this->setValueInputWidget(true);
		}
	}
}


void FEMConditionSetupDialog::on_disTypeBox_currentIndexChanged(int index)
{
	this->setValueInputWidget(index>0);
}

void FEMConditionSetupDialog::setValueInputWidget(bool is_button)
{
	if (is_button) // linear or direct
	{
		static_cast<QGridLayout*>(this->layout())->removeWidget(this->firstValueEdit);
		delete firstValueEdit;
		firstValueEdit = nullptr;
		directButton = new QPushButton("Calculate Values");
		connect(this->directButton, SIGNAL(pressed()), this, SLOT(directButton_pressed()));
		static_cast<QGridLayout*>(this->layout())->addWidget(directButton,5,1);
	}
	else	// constant or domain
	{
		static_cast<QGridLayout*>(this->layout())->removeWidget(this->directButton);
		delete directButton;
		directButton = nullptr;
		firstValueEdit = new QLineEdit("0");
		this->firstValueEdit->setValidator (_first_value_validator);
		static_cast<QGridLayout*>(this->layout())->addWidget(this->firstValueEdit,5,1);
	}
}

void FEMConditionSetupDialog::directButton_pressed()
{
	if (this->_mesh == nullptr)
	{
		const GeoLib::Polyline* line = dynamic_cast<const GeoLib::Polyline*>(_cond.getGeoObj());
		const std::vector<size_t> nodes = _cond.getDisNodes();
		const std::vector<double> values = _cond.getDisValues();
		LinearEditDialog dlg(*line, nodes, values);
		connect(&dlg, SIGNAL(transmitDisValues(std::vector< std::pair<std::size_t,double> >)),
				this, SLOT(addDisValues(std::vector< std::pair<std::size_t,double> >)));
		dlg.exec();
	}
	else
	{
		std::vector<MeshLib::Mesh*> msh_vec;
		msh_vec.push_back( const_cast<MeshLib::Mesh*>(this->_mesh) );
		CondFromRasterDialog dlg(msh_vec);
		//connect(&dlg, SIGNAL(directNodesWritten(std::string)), this, SLOT(direct_path_changed(std::string)));
		connect(&dlg, SIGNAL(transmitDisValues(std::vector< std::pair<std::size_t,double> >)),
				this, SLOT(addDisValues(std::vector< std::pair<std::size_t,double> >)));
		dlg.exec();
	}
}

void FEMConditionSetupDialog::addDisValues(std::vector< std::pair<size_t,double> > direct_values)
{
	_cond.setDisValues(direct_values);
	this->directButton->setText(QString::number(direct_values.size()) + " values added");
}

FEMCondition* FEMConditionSetupDialog::typeCast(const FEMCondition &cond)
{
	FEMCondition* new_cond(nullptr);
	switch(this->condTypeBox->currentIndex())
	{
		case 0:
			new_cond = new BoundaryCondition(cond);
			break;
		case 1:
			new_cond = new InitialCondition(cond);
			break;
		default:
			new_cond = new SourceTerm(cond);
	}
	return new_cond;
}

void FEMConditionSetupDialog::copyCondOnPoints()
{
	std::vector<FEMCondition*> conditions;
	if (_cond.getGeomType() == GeoLib::GEOTYPE::POLYLINE)
	{
		const GeoLib::Polyline* ply = dynamic_cast<const GeoLib::Polyline*>(_cond.getGeoObj());
		size_t nPoints = ply->getNumberOfPoints();
		for (size_t i=0; i<nPoints; i++)
		{
			FEMCondition* cond = new FEMCondition(_cond);
			cond->setGeoObj(nullptr);
			cond->setGeoType(GeoLib::GEOTYPE::POINT);
			cond->setGeoName(_cond.getAssociatedGeometryName() + "_Point" + std::to_string(ply->getPointID(i)));
			cond->clearDisValues();
			cond->setConstantDisValue((*ply->getPoint(i))[2]);
			conditions.push_back(this->typeCast(*cond));
		}
		emit createFEMCondition(conditions);
	}
	else if (_cond.getGeomType() == GeoLib::GEOTYPE::SURFACE)
	{
		const GeoLib::Surface* sfc = dynamic_cast<const GeoLib::Surface*>(_cond.getGeoObj());
		size_t nTriangles = sfc->getNTriangles();
		for (size_t i=0; i<nTriangles; i++)
		{
			const GeoLib::Triangle* tri = (*sfc)[i];
			for (size_t j=0; j<3; j++)
			{
				FEMCondition* cond = new FEMCondition(_cond);
				cond->setGeoObj(nullptr);
				cond->setGeoType(GeoLib::GEOTYPE::POINT);
				cond->setGeoName(_cond.getAssociatedGeometryName() + "_Point" + std::to_string((*tri)[j]));
				cond->clearDisValues();
				cond->setConstantDisValue((*tri->getPoint(j))[2]);
				conditions.push_back(this->typeCast(*cond));
			}
		}
		emit createFEMCondition(conditions);
	}
	else
		ERR("FEMConditionSetupDialog::copyCondOnPoints(): discerning GeoType.");
}

void FEMConditionSetupDialog::clearDisTypeBox()
{
	while (this->disTypeBox->count()>0)
		this->disTypeBox->removeItem(0);
}
