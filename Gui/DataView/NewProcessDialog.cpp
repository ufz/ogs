/**
 * \file NewProcessDialog.cpp
 * 2011/11/17 KR Initial implementation
 */

#include "NewProcessDialog.h"
#include "FEMEnums.h"
#include "ProcessInfo.h"


NewProcessDialog::NewProcessDialog(QDialog* parent)
: QDialog(parent)
{
	setupUi(this);
	setupDialog();
}

void NewProcessDialog::setupDialog()
{
	const std::list<std::string> process_names = FiniteElement::getAllProcessNames();
	for (std::list<std::string>::const_iterator it = process_names.begin(); it != process_names.end(); ++it)
		this->processTypeBox->addItem(QString::fromStdString(*it));

	const std::list<std::string> pv_names = FiniteElement::getAllPrimaryVariableNames();
	for (std::list<std::string>::const_iterator it = pv_names.begin(); it != pv_names.end(); ++it)
		this->pvTypeBox->addItem(QString::fromStdString(*it));
}

void NewProcessDialog::accept()
{
	ProcessInfo* info = new ProcessInfo();
	info->setProcessType(static_cast<FiniteElement::ProcessType>(this->processTypeBox->currentIndex() + 1));
	info->setProcessPrimaryVariable(static_cast<FiniteElement::PrimaryVariable>(this->pvTypeBox->currentIndex() + 1));

	emit addProcess(info);
	this->done(QDialog::Accepted);
}

void NewProcessDialog::reject()
{
	this->done(QDialog::Rejected);
}
