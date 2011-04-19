/**
 * \file RecentFiles.h
 * 5/11/2009 LB Initial implementation
 *
 */


#ifndef RECENTFILES_H
#define RECENTFILES_H

// ** INCLUDES **
#include <QObject>
#include <QMenu>
#include <QAction>

class QString;

/**
 * The RecentFiles class provides functionality to store informations about
 * recently used files (e.g. loaded or saved files).
 * Example Usage:
 * \code RecentFiles* recentFiles = new RecentFiles(this, SLOT(openRecentFile()), "recentFileList", "OpenGeoSys-5");
 * connect(this, SIGNAL(fileUsed(QString)), recentFiles, SLOT(setCurrentFile(QString)));
 * menu_File->addMenu( recentFiles->menu() ); \endcode
 * with:
 * \code void MainWindow::openRecentFile()
 * {
 *   QAction* action = qobject_cast<QAction*>(sender());
 *   if (action)
 *   loadFile(action->data().toString());
 * } \endcode
 */
class RecentFiles : public QObject
{
	Q_OBJECT

public:
	/**
	 * Constructor. Example Usage:
	 * \code RecentFiles recentFiles = new RecentFiles(this, SLOT(mySlot(QString)), "recentFileList"); \endcode
	 * \param parent The parent object. Normally the QMainWindow instance
	 * \param slot A slot on parent which is called when a recent file is clicked.
	 * Use this with Qts SLOT() macro!
	 * \param settingsName The setting key
	 * \param programName The name of the program. QSettings of one program
	 * should be stored with the same keys: QSettings("UFZ", programName)
	 */
	RecentFiles(QObject* parent, const char* slot, QString settingsName, QString programName);
	~RecentFiles();

	/// Returns the created menu. Add this menu to your QMainWindow menu.
	QMenu* menu();

public slots:
	/// Should be called from the application when a file was used.
	void setCurrentFile(const QString& filename);

private:
	/// Updates the recent files list and writes it to the settings.
	void updateRecentFileActions();

	/// Returns the filename from a full file path.
	QString strippedName(const QString& fullFileName);

	QMenu* _filesMenu;
	QString _currentFile;
	QString _settingsName;
	QString _programName;
	enum { _maxFiles = 5 };
	QAction* _fileActions[_maxFiles];

};

#endif // RECENTFILES_H
