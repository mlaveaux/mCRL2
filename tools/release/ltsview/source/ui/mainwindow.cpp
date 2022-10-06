// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <QApplication>
#include <QImageWriter>
#include <QSettings>
#include <QTimer>
#include <QSurfaceFormat>

#include "mcrl2/lts/lts_io.h"
#include "mainwindow.h"
#include "savepicturedialog.h"

MainWindow::MainWindow(QThread *atermThread):
  m_fileDialog("", this)
{
  m_ui.setupUi(this);

  m_ltsManager = new LtsManager(this, atermThread);
  m_markManager = new MarkManager(this, m_ltsManager);

  m_infoDock = new InfoDock(this, m_ltsManager, m_markManager);
  m_markDock = new MarkDock(this, m_markManager);
  m_simDock = new SimDock(this, m_ltsManager);
  m_settingsDock = new SettingsDock(this);
  m_settingsDialog = new SettingsDialog(this);

  
  QSurfaceFormat format = QSurfaceFormat(QSurfaceFormat::DebugContext);
  format.setMajorVersion(4);
  format.setMinorVersion(3);
  format.setProfile(QSurfaceFormat::CoreProfile);
  QSurfaceFormat::setDefaultFormat(format);

  m_glwidget = new GLWidget(nullptr, m_ltsManager, this);
  m_glwidget->setFormat(format);
  
  m_graphics_info_dialog = new GraphicsInfoDialog(this);
  setCentralWidget(m_glwidget);
  m_progressDialog = new QProgressDialog("", QString(), 0, 6, this);
  m_progressDialog->setMinimumDuration(0);

  m_ui.informationDock->setWidget(m_infoDock);
  m_ui.simulationDock->setWidget(m_simDock);
  m_ui.markDock->setWidget(m_markDock);
  m_ui.settingsDock->setWidget(m_settingsDock);

  m_ui.viewMenu->insertAction(m_ui.preferences, m_ui.informationDock->toggleViewAction());
  m_ui.viewMenu->insertAction(m_ui.preferences, m_ui.simulationDock->toggleViewAction());
  m_ui.viewMenu->insertAction(m_ui.preferences, m_ui.markDock->toggleViewAction());
  m_ui.viewMenu->insertAction(m_ui.preferences, m_ui.settingsDock->toggleViewAction());
  m_ui.viewMenu->insertSeparator(m_ui.preferences);

  mcrl2::log::logger::register_output_policy(m_logRelay);
  connect(&m_logRelay, SIGNAL(logMessage(QString, QDateTime, QString)), this, SLOT(logMessage(QString, QDateTime, QString)));

  connect(m_ui.open, SIGNAL(triggered()), this, SLOT(open()));
  connect(m_ui.openTrace, SIGNAL(triggered()), this, SLOT(openTrace()));
  connect(m_ui.exportBitmap, SIGNAL(triggered()), this, SLOT(exportBitmap()));
  connect(m_ui.exportText, SIGNAL(triggered()), this, SLOT(exportText()));
  connect(m_ui.exit, SIGNAL(triggered()), QApplication::instance(), SLOT(quit()));

  /// TODO: restore connection
  // connect(m_ui.resetViewpoint, SIGNAL(triggered()), m_ltsCanvas, SLOT(resetView()));
  connect(m_ui.zoomIntoAbove, SIGNAL(triggered()), m_ltsManager, SLOT(zoomInAbove()));
  connect(m_ui.zoomIntoBelow, SIGNAL(triggered()), m_ltsManager, SLOT(zoomInBelow()));
  connect(m_ui.zoomOut, SIGNAL(triggered()), m_ltsManager, SLOT(zoomOut()));

  connect(m_ui.displayStates, SIGNAL(triggered(bool)), &Settings::instance().displayStates, SLOT(setValue(bool)));
  connect(&Settings::instance().displayStates, SIGNAL(changed(bool)), m_ui.displayStates, SLOT(setChecked(bool)));
  connect(m_ui.displayTransitions, SIGNAL(triggered(bool)), &Settings::instance().displayTransitions, SLOT(setValue(bool)));
  connect(&Settings::instance().displayTransitions, SIGNAL(changed(bool)), m_ui.displayTransitions, SLOT(setChecked(bool)));
  connect(m_ui.displayBackpointers, SIGNAL(triggered(bool)), &Settings::instance().displayBackpointers, SLOT(setValue(bool)));
  connect(&Settings::instance().displayBackpointers, SIGNAL(changed(bool)), m_ui.displayBackpointers, SLOT(setChecked(bool)));
  connect(m_ui.displayWireframe, SIGNAL(triggered(bool)), &Settings::instance().displayWireframe, SLOT(setValue(bool)));
  connect(&Settings::instance().displayWireframe, SIGNAL(changed(bool)), m_ui.displayWireframe, SLOT(setChecked(bool)));

  connect(m_ui.preferences, SIGNAL(triggered()), m_settingsDialog, SLOT(show()));
  connect(m_ui.actionShowGraphicsInfo, SIGNAL(triggered()), m_graphics_info_dialog, SLOT(show()));

  connect(m_ltsManager, SIGNAL(loadingLts()), this, SLOT(loadingLts()));
  connect(m_ltsManager, SIGNAL(rankingStates()), this, SLOT(rankingStates()));
  connect(m_ltsManager, SIGNAL(clusteringStates()), this, SLOT(clusteringStates()));
  connect(m_ltsManager, SIGNAL(computingClusterInfo()), this, SLOT(computingClusterInfo()));
  connect(m_ltsManager, SIGNAL(positioningClusters()), this, SLOT(positioningClusters()));
  connect(m_ltsManager, SIGNAL(positioningStates()), this, SLOT(positioningStates()));
  connect(m_ltsManager, SIGNAL(ltsStructured()), this, SLOT(hideProgressDialog()));
  connect(m_ltsManager, SIGNAL(errorLoadingLts()), this, SLOT(hideProgressDialog()));
  connect(m_ltsManager, SIGNAL(startStructuring()), this, SLOT(startStructuring()));
  connect(m_ltsManager, SIGNAL(stopStructuring()), this, SLOT(stopStructuring()));

  connect(m_ltsManager, SIGNAL(selectionChanged()), this, SLOT(selectionChanged()));
  connect(m_ltsManager, SIGNAL(ltsZoomed(LTS *)), this, SLOT(zoomChanged()));
  m_ui.zoomIntoAbove->setEnabled(false);
  m_ui.zoomIntoBelow->setEnabled(false);
  m_ui.zoomOut->setEnabled(false);


  /// TODO: Restore connections
  // connect(m_ltsCanvas, SIGNAL(renderingStarted()), this, SLOT(startRendering()));
  // connect(m_ltsCanvas, SIGNAL(renderingFinished()), this, SLOT(clearStatusBar()));
  // connect(m_ltsCanvas, SIGNAL(renderingFinished()), this, SLOT(updateGraphicsInfo()));

  QSettings settings("mCRL2", "LTSView");
  restoreGeometry(settings.value("geometry").toByteArray());
  restoreState(settings.value("windowState").toByteArray());
}

MainWindow::~MainWindow()
{
  mcrl2::log::logger::unregister_output_policy(m_logRelay);
}

void MainWindow::updateGraphicsInfo(){ 
  // if (!m_ltsCanvas) return;
  // m_graphics_info_dialog->lbl_info->setText(m_ltsCanvas->graphics_info);
}

void MainWindow::closeEvent(QCloseEvent * /*event*/)
{
  {
    QSettings settings("mCRL2", "LTSView");
    settings.setValue("geometry", saveGeometry());
    settings.setValue("windowState", saveState());
  }

  // It's possible we're in the middle of an LTS restructuring operation that can't be stopped in any friendly way, so hard-terminate the program.
  _exit(0);
}

void MainWindow::open(QString filename)
{
  loadingLts();
  m_progressDialog->setWindowTitle("Opening file");

  if (m_ltsManager->openLts(filename))
  {
    setWindowTitle(filename + " - LTSView");
  }
}

void MainWindow::open()
{
  QString filter = QString("All supported files (") + QString::fromStdString(mcrl2::lts::detail::lts_extensions_as_string(" ")) + ");;All files (*.*)";
  QString filename = m_fileDialog.getOpenFileName("Open LTS", filter);
  if (filename.isNull())
  {
    return;
  }

  open(filename);
}

void MainWindow::openTrace()
{
  if (!m_ltsManager->lts())
  {
    return;
  }

  QString filename = m_fileDialog.getOpenFileName("Open Trace", "Traces (*.trc);;All files (*.*)");
  if (filename.isNull())
  {
    return;
  }

  m_ltsManager->loadTrace(filename);
}

void MainWindow::exportBitmap()
{
  QString filter = "Images (";
  QList<QByteArray> formats = QImageWriter::supportedImageFormats();
  for(QList<QByteArray>::iterator i = formats.begin(); i != formats.end(); i++)
  {
    if(i != formats.begin())
    {
      filter += " ";
    }
    filter += "*." + *i;
  }
  filter += ")";

  QString filename = m_fileDialog.getSaveFileName("Save picture", filter);
  if (filename.isNull())
  {
    return;
  }
/// TODO: Restore connection
  // SavePictureDialog dialog(this, m_ltsCanvas, filename);
  // connect(&dialog, SIGNAL(statusMessage(QString)), this, SLOT(setStatusBar(QString)));
  // dialog.exec();
  clearStatusBar();
}

void MainWindow::exportText()
{
  m_glwidget->setRoot(m_ltsManager->lts()->getInitialState()->getCluster());
  QString filename = m_fileDialog.getSaveFileName("Save text", "All files (*.*)");
  if (filename.isNull())
  {
    return;
  }
  /// TODO: Restore functionality
  // m_ltsCanvas->exportToText(filename);
}

void MainWindow::setProgress(int phase, QString message)
{
  mCRL2log(mcrl2::log::verbose) << "LTSView MainWindow::setProgress(" << phase << ", " << message.toStdString() << ")" << std::endl;
  if (!m_progressDialog->isVisible())
  {
    m_progressDialog->setWindowTitle("Structuring LTS");
    m_progressDialog->show();
  }
  m_progressDialog->setLabelText(message);
  m_progressDialog->setValue(phase);
}

void MainWindow::selectionChanged()
{
  m_ui.zoomIntoAbove->setEnabled(m_ltsManager->selectedCluster() != 0);
  m_ui.zoomIntoBelow->setEnabled(m_ltsManager->selectedCluster() != 0);
}
