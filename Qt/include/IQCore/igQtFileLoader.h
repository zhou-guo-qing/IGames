/**
 * @class   igQtFileLoader
 * @brief   igQtFileLoader's brief
 */

#pragma once
#include "iGameSceneManager.h"

#include <QFileDialog>
#include <QString>
#include <IQCore/igQtExportModule.h>

using namespace iGame;

class IG_QT_MODULE_EXPORT igQtFileLoader : public QObject
{
	Q_OBJECT
public:
	igQtFileLoader(QObject* parent = nullptr);
	~igQtFileLoader() override;

public:
	void LoadFile();
	void OpenFile(const std::string& fileName);
	void SaveFile(); //保存
	void SaveFileAs();//另存为
	void SaveCurrentFileToRecentFile(QString file_name);
	void AddCurrentFileToRecentFilePath(QString lastPath);
	void InitRecentFilePaths();
	void InitRecentFileActions(std::vector<QString>);
	void UpdateRecentActionList();
	void UpdateIniFileInfo();
	QList<QAction*> GetRecentActionList() { return this->recentFileActionList; };

signals:
	void NewModel(DataObject::Pointer obj, ItemSource source);

	void FinishReading();
	void EmitMakeCurrent();
	void EmitDoneCurrent();

	void AddFileToModelList(QString file_name);

	void LoadAnimationFile(std::vector<float>& timeValues);

protected:
	QList<QAction*> recentFileActionList;
	int maxFileNr = 10;
    SceneManager::Pointer m_SceneManager;
};
