/**
 * \file SHPImportDialog.h
 * 25/01/2010 KR Initial implementation
 */

#ifndef SHPIMPORTDIALOG_H
#define SHPIMPORTDIALOG_H

#include <QDialog>
#include <QtGui/QMainWindow>

class SHPInterface;
class GEOModels;

class QDialogButtonBox;
class QFileInfo;
class QGridLayout;
class QLabel;
class QLineEdit;
class QRadioButton;
class QVBoxLayout;

/**
 * \brief Dialog for selecting which information should be loaded from a shape file.
 */
class SHPImportDialog : public QDialog
{
	Q_OBJECT

public:
	/// Constructor
	SHPImportDialog(std::string filename, GEOModels* geoModels, QDialog* parent = 0);
	~SHPImportDialog();

	QDialogButtonBox* _buttonBox; /// The buttons used in this dialog.

private:
	/// Constructs a dialog window based on the information found in the selected shape file
	void setupDialog();

	QGridLayout* _layout;
	QLabel* _shpContentLabel;
	QLabel* _nameLabel;
	QLineEdit* _listName;
	QRadioButton* _choice1, * _choice2;
	std::string _filename;
	short _fileType;
	SHPInterface* _shpInterface;

private slots:
	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

signals:
	void shpLoaded(QString);
};

#endif //SHPIMPORTDIALOG_H
