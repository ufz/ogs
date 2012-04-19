/**
 * \file FEMConditionSetupDialog.cpp
 * 2011/11/07 KR Initial implementation
 */

#include "FEMConditionSetupDialog.h"
#include "OGSError.h"
#include "FEMEnums.h"
#include "ProjectData.h"
#include "StrictDoubleValidator.h"
#include "StringTools.h"
#include "LinearEditDialog.h"
#include "CondFromRasterDialog.h"

#include "BoundaryCondition.h"
#include "InitialCondition.h"
#include "SourceTerm.h"

FEMConditionSetupDialog::FEMConditionSetupDialog(const std::string &associated_geometry,
												 const GEOLIB::GEOTYPE type,
												 const std::string &geo_name,
												 const GEOLIB::GeoObject* const geo_object,
												 bool  on_points,
												 QDialog* parent)
: QDialog(parent), _cond(associated_geometry, FEMCondition::UNSPECIFIED), _set_on_points(on_points), 
  _combobox(NULL), directButton(NULL), _mesh(NULL), _first_value_validator(NULL)
{
	_cond.setGeoType(type);
	_cond.setGeoName(geo_name);
	_cond.setGeoObj(geo_object);

	setupUi(this);
	setupDialog();
}

FEMConditionSetupDialog::FEMConditionSetupDialog(const FEMCondition &cond, QDialog* parent)
	: QDialog(parent), _cond(cond), _set_on_points(false), _combobox(NULL), directButton(NULL), 
	_mesh(NULL), _first_value_validator(NULL)
{
	setupUi(this);
	setupDialog();
	setValuesFromCond();
}

FEMConditionSetupDialog::FEMConditionSetupDialog(const std::string &name, const MeshLib::CFEMesh* mesh, QDialog* parent)
: QDialog(parent), _cond(name, FEMCondition::UNSPECIFIED), _set_on_points(false),  _combobox(NULL), directButton(NULL),
  _mesh(mesh), _first_value_validator(NULL)
{
	_cond.setGeoType(GEOLIB::INVALID);
	_cond.setGeoName(name);
	_cond.setGeoObj(NULL);

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
	if (_cond.getGeoType() != GEOLIB::INVALID)
	{
		this->disTypeBox->addItem("Constant (Direchlet)");
		if (_cond.getGeoType() == GEOLIB::POLYLINE)
			this->disTypeBox->addItem("Linear (Direchlet)");

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
		static_cast<QGridLayout*>(this->layout())->removeWidget(this->firstValueEdit);
		directButton = new QPushButton("Calculate Values");
		static_cast<QGridLayout*>(this->layout())->addWidget(directButton,5,1);
		connect(this->directButton, SIGNAL(pressed()), this, SLOT(directButton_pressed()));
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

	int idx(0);
	if (_cond.getCondType() == FEMCondition::INITIAL_CONDITION)
		this->condTypeBox->setCurrentIndex(1);
	else if (_cond.getCondType() == FEMCondition::SOURCE_TERM)
	{
		this->condTypeBox->setCurrentIndex(2);
		on_condTypeBox_currentIndexChanged(2);
	}

	if (_cond.getGeoType() != GEOLIB::INVALID)
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
		size_t nValues (_cond.getDisValues().size());
		this->directButton->setText(QString::number(static_cast<int>(_cond.getDisValues().size())) + " values");
	}
}


void FEMConditionSetupDialog::accept()
{
	_cond.setProcessType(static_cast<FiniteElement::ProcessType>(this->processTypeBox->currentIndex() + 1));
	_cond.setProcessPrimaryVariable(static_cast<FiniteElement::PrimaryVariable>(this->pvTypeBox->currentIndex() + 1));

	if (_cond.getGeoType() != GEOLIB::INVALID)
	{
		if (condTypeBox->currentIndex()>1)
		{
			if (this->disTypeBox->currentIndex()>0) 
				_cond.setProcessDistributionType(FiniteElement::LINEAR_NEUMANN);
			else 
			{
				_cond.setProcessDistributionType(FiniteElement::CONSTANT_NEUMANN);
				_cond.setConstantDisValue(this->firstValueEdit->text().toDouble());
			}
		}
		else
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
		_cond.setProcessDistributionType(FiniteElement::DIRECT);
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
	if (_cond.getGeoType() != GEOLIB::INVALID)
	{
		if (index>1) // source terms selected
		{
			while (this->disTypeBox->count()>0)
				this->disTypeBox->removeItem(0);
			this->disTypeBox->addItem("Constant (Neumann)");
			if (_cond.getGeoType() == GEOLIB::POLYLINE)
				this->disTypeBox->addItem("Linear (Neumann)");
		}
		else
		{
			while (this->disTypeBox->count()>0)
				this->disTypeBox->removeItem(0);
			this->disTypeBox->addItem("Constant (Direchlet)");
			if (_cond.getGeoType() == GEOLIB::POLYLINE)
				this->disTypeBox->addItem("Linear (Direchlet)");
		}
	}
}


void FEMConditionSetupDialog::on_disTypeBox_currentIndexChanged(int index)
{
	if (index>0) // linear
	{
		static_cast<QGridLayout*>(this->layout())->removeWidget(this->firstValueEdit);
		directButton = new QPushButton("Calculate Values");
		connect(this->directButton, SIGNAL(pressed()), this, SLOT(directButton_pressed()));
		static_cast<QGridLayout*>(this->layout())->addWidget(directButton,5,1);
	}
	else	// constant
	{
		static_cast<QGridLayout*>(this->layout())->removeWidget(this->directButton);
		delete directButton;
		directButton = NULL;
		static_cast<QGridLayout*>(this->layout())->addWidget(this->firstValueEdit,5,1);
	}

}

void FEMConditionSetupDialog::directButton_pressed()
{
	if (this->_mesh == NULL)
	{
		const GEOLIB::Polyline* line = dynamic_cast<const GEOLIB::Polyline*>(_cond.getGeoObj());
		const std::vector<size_t> nodes = _cond.getDisNodes();
		const std::vector<double> values = _cond.getDisValues();
		LinearEditDialog dlg(*line, nodes, values);
		connect(&dlg, SIGNAL(transmitDisValues(std::vector< std::pair<size_t,double> >)), 
				this, SLOT(addDisValues(std::vector< std::pair<size_t,double> >)));
		dlg.exec();
	}
	else
	{
		std::map<std::string, MeshLib::CFEMesh*> msh_map;
		msh_map.insert( std::pair<std::string, MeshLib::CFEMesh*>(this->_cond.getGeoName(), const_cast<MeshLib::CFEMesh*>(this->_mesh)) );
		CondFromRasterDialog dlg(msh_map);
		//connect(&dlg, SIGNAL(directNodesWritten(std::string)), this, SLOT(direct_path_changed(std::string)));
		connect(&dlg, SIGNAL(transmitDisValues(std::vector< std::pair<size_t,double> >)), 
				this, SLOT(addDisValues(std::vector< std::pair<size_t,double> >)));
		dlg.exec();
	}
}

void FEMConditionSetupDialog::addDisValues(std::vector< std::pair<size_t,double> > direct_values)
{
	_cond.setDisValues(direct_values);
	this->directButton->setText(QString::number(direct_values.size()) + " values added");
	//this->directButton->setEnabled(false);
}

FEMCondition* FEMConditionSetupDialog::typeCast(const FEMCondition &cond)
{
	FEMCondition* new_cond(NULL);
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
	if (_cond.getGeoType() == GEOLIB::POLYLINE)
	{
		const GEOLIB::Polyline* ply = dynamic_cast<const GEOLIB::Polyline*>(_cond.getGeoObj());
		size_t nPoints = ply->getNumberOfPoints();
		for (size_t i=0; i<nPoints; i++)
		{
			FEMCondition* cond = new FEMCondition(_cond);
			cond->setGeoObj(NULL);
			cond->setGeoType(GEOLIB::POINT);
			cond->setGeoName(_cond.getAssociatedGeometryName() + "_Point" + number2str(ply->getPointID(i)));
			cond->clearDisValues();
			cond->setConstantDisValue((*ply->getPoint(i))[2]);
			conditions.push_back(this->typeCast(*cond));
		}
		emit createFEMCondition(conditions);
	}
	else if (_cond.getGeoType() == GEOLIB::SURFACE)
	{
		const GEOLIB::Surface* sfc = dynamic_cast<const GEOLIB::Surface*>(_cond.getGeoObj());
		size_t nTriangles = sfc->getNTriangles();
		for (size_t i=0; i<nTriangles; i++)
		{
			const GEOLIB::Triangle* tri = (*sfc)[i];
			for (size_t j=0; j<3; j++)
			{
				FEMCondition* cond = new FEMCondition(_cond);
				cond->setGeoObj(NULL);
				cond->setGeoType(GEOLIB::POINT);
				cond->setGeoName(_cond.getAssociatedGeometryName() + "_Point" + number2str((*tri)[j]));
				cond->clearDisValues();
				cond->setConstantDisValue((*tri->getPoint(j))[2]);
				conditions.push_back(this->typeCast(*cond));
			}
		}	
		emit createFEMCondition(conditions);
	}
	else
		std::cout << "Error discerning GeoType ..." << std::endl;
}
