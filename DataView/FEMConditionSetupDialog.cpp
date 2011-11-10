/**
 * \file FEMConditionSetupDialog.cpp
 * 2011/11/07 KR Initial implementation
 */

#include "FEMConditionSetupDialog.h"
#include "OGSError.h"
#include "FEMEnums.h"
#include "ProjectData.h"
#include "StrictDoubleValidator.h"

#include "BoundaryCondition.h"
#include "InitialCondition.h"
#include "SourceTerm.h"

FEMConditionSetupDialog::FEMConditionSetupDialog(const std::string &associated_geometry, 
												 const GEOLIB::GEOTYPE type, 
												 const std::string &geo_name,
												 const GEOLIB::GeoObject* const geo_object, 
												 QDialog* parent)
: QDialog(parent), _cond(associated_geometry, FEMCondition::UNSPECIFIED), _secondValueEdit(NULL),
  _first_value_validator(NULL), _second_value_validator(NULL)
{
	_cond.setGeoType(type);
	_cond.setGeoName(geo_name);
	_cond.setGeoObj(geo_object);

	setupUi(this);

	setupDialog();


}

FEMConditionSetupDialog::FEMConditionSetupDialog(FEMCondition &cond, QDialog* parent)
	: QDialog(parent), _cond(cond), _secondValueEdit(NULL), 
	  _first_value_validator(NULL), _second_value_validator(NULL)
{
	setupDialog();
}

FEMConditionSetupDialog::~FEMConditionSetupDialog()
{
	delete _secondValueEdit;
	delete _first_value_validator;
	delete _second_value_validator;
}

void FEMConditionSetupDialog::setupDialog()
{
	if (_cond.getGeoType() == GEOLIB::POLYLINE)
	{
		this->disTypeBox->addItem("Linear (Direchlet)");
		this->disTypeBox->addItem("Linear (Neumann)");
	}
	_first_value_validator = new StrictDoubleValidator(-1e+10, 1e+10, 5);
	_second_value_validator = new StrictDoubleValidator(-1e+10, 1e+10, 5);
	this->firstValueEdit->setText("0");
	this->firstValueEdit->setValidator (_first_value_validator);

	const std::list<std::string> process_names = getAllProcessNames();
	for (std::list<std::string>::const_iterator it = process_names.begin(); it != process_names.end(); ++it)
		this->processTypeBox->addItem(QString::fromStdString(*it));

	const std::list<std::string> pv_names = getAllPrimaryVariableNames();
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
	_cond.setProcessType(static_cast<ProcessType>(this->processTypeBox->currentIndex() + 1));
	_cond.setProcessPrimaryVariable(static_cast<PrimaryVariable>(this->pvTypeBox->currentIndex() + 1));

	switch(this->disTypeBox->currentIndex())
	{
		case 0:
			_cond.setProcessDistributionType(FiniteElement::CONSTANT);
			break;
		case 1:
			_cond.setProcessDistributionType(FiniteElement::CONSTANT_NEUMANN);
			break;
		case 2:
			_cond.setProcessDistributionType(FiniteElement::LINEAR);
			break;
		case 3:
			_cond.setProcessDistributionType(FiniteElement::LINEAR_NEUMANN);
			break;
		default:
			_cond.setProcessDistributionType(FiniteElement::INVALID_DIS_TYPE);
	}
	
	std::vector<double> dis_values;
	dis_values.push_back(strtod(this->firstValueEdit->text().toStdString().c_str(), 0));
	if (this->_secondValueEdit)
		dis_values.push_back(strtod(this->_secondValueEdit->text().toStdString().c_str(), 0));
	_cond.setDisValue(dis_values);

	FEMCondition* new_cond(NULL);
	switch(this->condTypeBox->currentIndex())
	{
		case 0:
			new_cond = new BoundaryCondition(_cond);
			break;
		case 1:
			new_cond = new InitialCondition(_cond);
			break;
		default:
			new_cond = new SourceTerm(_cond);
	}

	emit addFEMCondition(new_cond);
	this->done(QDialog::Accepted);
}

void FEMConditionSetupDialog::reject()
{
	this->done(QDialog::Rejected);
}

/*
void FEMConditionSetupDialog::on_condTypeBox_currentIndexChanged(int index)
{
	if (index==1)
		this->geoNameBox->addItem("Domain");
	// remove "Domain" if IC is unselected
}
*/

void FEMConditionSetupDialog::on_disTypeBox_currentIndexChanged(int index)
{
	if (index>1) // linear
	{
		if (!_secondValueEdit)
		{
			_secondValueEdit = new QLineEdit("0");
			_secondValueEdit->setValidator(_second_value_validator);
			static_cast<QGridLayout*>(this->layout())->addWidget(_secondValueEdit,5,1) ;
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