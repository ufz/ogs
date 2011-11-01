/**
 * \file SetNameDialog.cpp
 * 2011/10/26 KR Initial implementation
 */

#include <SetNameDialog.h>
#include <QDialogButtonBox>
#include <QDialogButtonBox>
#include <QLabel>
#include <QLineEdit>
#include <QVBoxLayout>

SetNameDialog::SetNameDialog(const std::string &parent_name, const std::string &object_type_name, size_t id, const std::string &old_name = "", QDialog* parent) :
	_parent_name(parent_name), _object_type_name(object_type_name), _id(id)
{
	setupDialog(old_name);
	show();
}

SetNameDialog::~SetNameDialog()
{
	delete _buttonBox;
	delete _layout;
	delete _new_name;
	delete _txt_label;
}

void SetNameDialog::setupDialog(const std::string &old_name)
{
	_layout = new QVBoxLayout(this);
	QString dialog_text("Please enter a name for " + QString::fromStdString(_object_type_name) + " #" + QString::number(_id));
	_txt_label = new QLabel(this);
	_txt_label->setText(dialog_text);
	_new_name = new QLineEdit(QString::fromStdString(old_name));

	setWindowTitle("Set name...");
	_layout->addWidget( _txt_label );
	_layout->addWidget( _new_name );
	_buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
	connect(_buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
	connect(_buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
	_layout->addWidget( _buttonBox );

	setLayout(_layout);
}

void SetNameDialog::accept()
{
	emit requestNameChange(_parent_name, _object_type_name, _id, _new_name->text().toStdString());	
	this->done(QDialog::Accepted);
}

void SetNameDialog::reject()
{
	this->done(QDialog::Rejected);
}
