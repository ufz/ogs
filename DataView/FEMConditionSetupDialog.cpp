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
: QDialog(parent), _cond(associated_geometry, FEMCondition::UNSPECIFIED), _set_on_points(on_points), _secondValueEdit(NULL),
  directButton(NULL), _mesh(NULL), _first_value_validator(NULL), _second_value_validator(NULL)
{
	_cond.setGeoType(type);
	_cond.setGeoName(geo_name);
	_cond.setGeoObj(geo_object);

	setupUi(this);

	setupDialog();
}

FEMConditionSetupDialog::FEMConditionSetupDialog(FEMCondition &cond, QDialog* parent)
	: QDialog(parent), _cond(cond), _secondValueEdit(NULL), directButton(NULL),
	  _first_value_validator(NULL), _second_value_validator(NULL)
{
	setupDialog();
}

FEMConditionSetupDialog::FEMConditionSetupDialog(const std::string &name, const MeshLib::CFEMesh* mesh, QDialog* parent)
: QDialog(parent), _cond(name, FEMCondition::UNSPECIFIED), _set_on_points(false), _secondValueEdit(NULL),
  directButton(NULL), _mesh(mesh), _first_value_validator(NULL), _second_value_validator(NULL)
{
	_cond.setGeoType(GEOLIB::INVALID);
	_cond.setGeoName(name);
	_cond.setGeoObj(NULL);

	setupUi(this);

	setupDialog();
}

FEMConditionSetupDialog::~FEMConditionSetupDialog()
{
	delete directButton;
	delete _secondValueEdit;
	delete _first_value_validator;
	delete _second_value_validator;
}

void FEMConditionSetupDialog::setupDialog()
{
	if (_cond.getGeoType() != GEOLIB::INVALID)
	{
		this->disTypeBox->addItem("Constant (Direchlet)");
		if (_cond.getGeoType() == GEOLIB::POLYLINE)
			this->disTypeBox->addItem("Linear (Direchlet)");

		_first_value_validator = new StrictDoubleValidator(-1e+10, 1e+10, 5);
		_second_value_validator = new StrictDoubleValidator(-1e+10, 1e+10, 5);
		this->firstValueEdit->setText("0");
		this->firstValueEdit->setValidator (_first_value_validator);
	}
	else	// direct on mesh
	{
		this->valueLabel->setText("DirectNode File");
		directButton = new QPushButton("Calculate Values");
		static_cast<QGridLayout*>(this->layout())->addWidget(directButton,6,1) ;
		connect(this->directButton, SIGNAL(pressed()), this, SLOT(directButton_pressed()));
		this->disTypeBox->addItem("Direct");
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

void FEMConditionSetupDialog::accept()
{
	_cond.setProcessType(static_cast<FiniteElement::ProcessType>(this->processTypeBox->currentIndex() + 1));
	_cond.setProcessPrimaryVariable(static_cast<FiniteElement::PrimaryVariable>(this->pvTypeBox->currentIndex() + 1));

	if (_cond.getGeoType() != GEOLIB::INVALID)
	{
		//QString dis_type_text = this->disTypeBox->currentText();
		if (condTypeBox->currentIndex()>1)
		{
			if (this->disTypeBox->currentIndex()>0) 
				_cond.setProcessDistributionType(FiniteElement::LINEAR_NEUMANN);
			else 
				_cond.setProcessDistributionType(FiniteElement::CONSTANT_NEUMANN);
		}
		else
		{
			if (this->disTypeBox->currentIndex()>0) 
				_cond.setProcessDistributionType(FiniteElement::LINEAR);
			else 
				_cond.setProcessDistributionType(FiniteElement::CONSTANT);
		}

		std::vector<double> dis_values;
		dis_values.push_back(strtod(this->firstValueEdit->text().toStdString().c_str(), 0));
		if (this->_secondValueEdit)
			dis_values.push_back(strtod(this->_secondValueEdit->text().toStdString().c_str(), 0));
		_cond.setDisValues(dis_values);
	}
	else	// direct on mesh
	{
		_cond.setProcessDistributionType(FiniteElement::DIRECT);
		std::string direct_node_path = this->firstValueEdit->text().toStdString();
		_cond.setDirectFileName(direct_node_path);
		std::vector< std::pair<size_t, double> > node_values;
		SourceTerm::getDirectNodeValues(direct_node_path, node_values);
		_cond.setDisValues(node_values);
	}

	if (!_set_on_points)
		emit addFEMCondition(this->typeCast(_cond));
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
		if (!_secondValueEdit)
		{
			_secondValueEdit = new QLineEdit("0");
			_secondValueEdit->setValidator(_second_value_validator);
			static_cast<QGridLayout*>(this->layout())->addWidget(_secondValueEdit,6,1) ;
		}
	}
	else	// constant
	{
		if (_secondValueEdit)
		{
			static_cast<QGridLayout*>(this->layout())->removeWidget(_secondValueEdit);
			delete _secondValueEdit;
			_secondValueEdit = NULL;
		}
	}

}

void FEMConditionSetupDialog::directButton_pressed()
{
	std::map<std::string, MeshLib::CFEMesh*> msh_map;
	msh_map.insert( std::pair<std::string, MeshLib::CFEMesh*>(this->_cond.getGeoName(), const_cast<MeshLib::CFEMesh*>(this->_mesh)) );
	CondFromRasterDialog dlg(msh_map);
	connect(&dlg, SIGNAL(directNodesWritten(std::string)), this, SLOT(direct_path_changed(std::string)));
	dlg.exec();
}

void FEMConditionSetupDialog::direct_path_changed(std::string path)
{
	this->firstValueEdit->setText(QString::fromStdString(path));
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
			cond->setDisValue((*ply->getPoint(i))[2]);
			emit addFEMCondition(this->typeCast(*cond));
		}
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
				cond->setDisValue((*tri->getPoint(j))[2]);
				emit addFEMCondition(this->typeCast(*cond));
			}
		}	
	}
	else
		std::cout << "Error discerning GeoType ..." << std::endl;
}
