// Author(s): Ruben Vink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <QtOpenGL>
#include "glutil.h"
#include "mcrl2/utilities/logger.h"

#include <QWidget>

#include "glwidget.h"
#include "arcballcamera.h"
#include "glscenegraph.h"

void GLWidget::determineActiveTool(QMouseEvent* event, bool useModifiers)
{
  if (useModifiers && (event->modifiers() & Qt::ControlModifier) &&
      (event->modifiers() & Qt::ShiftModifier))
  {
    setActiveTool(ZoomTool);
  }
  else if (useModifiers && event->modifiers() & Qt::ControlModifier)
  {
    setActiveTool(PanTool);
  }
  else if ((event->buttons() & Qt::MidButton) ||
           ((event->buttons() & Qt::LeftButton) &&
            (event->buttons() & Qt::RightButton)))
  {
    setActiveTool(ZoomTool);
  }
  else if ((useModifiers && event->modifiers() & Qt::ShiftModifier) ||
           event->buttons() & Qt::RightButton)
  {
    setActiveTool(RotateTool);
  }
  else
  {
    setActiveTool(SelectTool);
  }
}

void GLWidget::setActiveTool(Tool tool)
{
  m_activeTool = tool;
  if (tool == PanTool)
  {
    setCursor(m_panCursor);
  }
  else if (tool == ZoomTool)
  {
    setCursor(m_zoomCursor);
  }
  else if (tool == RotateTool)
  {
    setCursor(m_rotateCursor);
  }
  else
  {
    setCursor(m_selectCursor);
  }
}

/// TODO: Use AABBs from the vistree to query hit object
GLWidget::Selection GLWidget::selectObject(QPoint position)
{
  return Selection();
}


GLWidget::GLWidget(Cluster* root, LtsManager* ltsManager, QWidget* parent)
    : QOpenGLWidget(parent), m_ltsManager(ltsManager), m_scene(root)
{
  connect(&Settings::instance().stateRankStyleCyclic, SIGNAL(changed(bool)),
          this, SLOT(onClusterStyleChanged()));
  connect(&Settings::instance().fsmStyle, SIGNAL(changed(bool)), this,
          SLOT(onClusterStyleChanged()));
  connect(&Settings::instance().statePosStyleMultiPass, SIGNAL(changed(bool)),
          this, SLOT(onClusterStyleChanged()));
  connect(&Settings::instance().clusterVisStyleTubes, SIGNAL(changed(bool)),
          this, SLOT(onClusterStyleChanged()));
}

void GLWidget::setRoot(Cluster* root)
{
  mCRL2log(mcrl2::log::debug) << "setRoot" << std::endl;
  m_scene.m_clusterRoot = root;
  if (m_initialised && QOpenGLContext::currentContext())
    m_scene.rebuildScene();
  //if (m_initialised) makeCurrent();
}

GLWidget::~GLWidget(){};

/**
 * @brief Initialises the OpenGL context.
 */
void GLWidget::initializeGL()
{
  // Check whether context creation succeeded and print the OpenGL major.minor
  // version.
  if (isValid())
  {

    // Check the minimum run-time requirement; the pair ordering is
    // lexicographical.
    QPair<int, int> version = format().version();
#ifndef __APPLE__
    QPair<int, int> required(3, 3);
    if (version < required)
    {
      // Print a message to the console and show a message box.
      std::stringstream message;

      message << "The runtime version of OpenGL (" << version.first << "."
              << version.second
              << ") is below the least supported version of OpenGL ("
              << required.first << "." << required.second << ").";
      mCRL2log(mcrl2::log::error) << message.str().c_str() << "\n";

      QMessageBox box(QMessageBox::Warning, "Unsupported OpenGL Version",
                      message.str().c_str(), QMessageBox::Ok);
      box.exec();

      throw mcrl2::runtime_error("Unsupported OpenGL version.");
    }
    else
#endif
    {

      QOpenGLContext::currentContext()->setFormat(
          QSurfaceFormat::defaultFormat());
      glEnable(GL_MULTISAMPLE);
      glCheckError();
      glEnable(GL_BLEND);
      glCheckError();
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glCheckError();

      mCRL2log(mcrl2::log::verbose) << "Created an OpenGL " << version.first
                                    << "." << version.second << " context.\n";
      //makeCurrent();
    }
  }
  else
  {
    mcrl2::runtime_error("Context creation failed.");
  }

  // Enable real-time logging of OpenGL errors when the GL_KHR_debug extension
  // is available. Ruben: Disabled because this makes the UI unusable with -d
  // flag
#if ENABLE_OPENGL_DEBUG_LOG
  m_logger = new QOpenGLDebugLogger(this);
  if (m_logger->initialize())
  {
    connect(m_logger, &QOpenGLDebugLogger::messageLogged, this,
            &GLWidget::logMessage);
    m_logger->startLogging();
  }
  else
  {
    mCRL2log(mcrl2::log::debug) << "QOpenGLDebugLogger initialisation failed\n";
  }
#else
  mCRL2log(mcrl2::log::debug)
      << "QOpenGLDebugLogger is disabled manually in \"glwidget.cpp\". Change "
         "\"#define ENABLE_OPENGL_DEBUG_LOG false\" to \"true\" if needed.\n";
#endif

  QOpenGLContext* ctx = QOpenGLContext::currentContext();
  mCRL2log(mcrl2::log::debug)
      << "GLWidget initializeGL: ctx: " << ctx << std::endl;
  m_scene.initializeScene();
  m_initialised = true;
}

void GLWidget::logMessage(const QOpenGLDebugMessage& debugMessage)
{
  mCRL2log(mcrl2::log::debug)
      << "OpenGL: " << debugMessage.message().toStdString() << "\n";
  std::cout << "OpenGL logger: " << debugMessage.message().toStdString()
            << "\n";
}

void GLWidget::paintGL()
{
  QPainter* painter = new QPainter(this);
  painter->beginNativePainting();
  QOpenGLFunctions* f = QOpenGLContext::currentContext()->functions();
  // Cull polygons that are facing away (back) from the camera, where their
  // front is defined as counter clockwise by default, see glFrontFace, meaning
  // that the vertices that make up a triangle should be oriented counter
  // clockwise to show the triangle.
  f->glFrontFace(GL_CW);

  // Enable depth testing, so that we don't have to care too much about
  // rendering in the right order. 
  // f->glEnable(GL_DEPTH_TEST);

  // Change the alpha blending function to make an alpha of 1 opaque and 0 fully
  // transparent. f->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  // f->glEnable(GL_BLEND);

  // Enable multisample antialiasing.
  f->glEnable(GL_MULTISAMPLE);

  f->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  f->glClearColor(Settings::instance().backgroundColor.value().redF(),
                  Settings::instance().backgroundColor.value().greenF(),
                  Settings::instance().backgroundColor.value().blueF(), 1.0);
  m_scene.renderScene();
  painter->endNativePainting();
};

void GLWidget::resizeGL(int width, int height)
{
  QOpenGLFunctions* f = QOpenGLContext::currentContext()->functions();
  f->glViewport(0, 0, width, height);
  m_width = width; m_height = height;
  m_scene.resizeScene(width, height);
};

void GLWidget::mousePressEvent(QMouseEvent* e)
{
  mCRL2log(mcrl2::log::debug) << "CLIKCKK" << std::endl;
  determineActiveTool(e, true);

  if (m_activeTool == SelectTool)
  {
    if (m_ltsManager->lts())
    {
      Selection selection = selectObject(e->pos());
      if (selection.state)
      {
        m_ltsManager->selectState(selection.state);
      }
      else if (selection.cluster)
      {
        m_ltsManager->selectCluster(selection.cluster);
      }
      else
      {
        m_ltsManager->unselect();
      }
    }
  }
  else
  {
    m_dragging = true;
    m_lastMousePosition = e->pos();
  }
  e->accept();
};

void GLWidget::mouseReleaseEvent(QMouseEvent* e){
  determineActiveTool(e, false);
  m_dragging = m_activeTool != SelectTool;
  e->accept();
  repaint();
};

void GLWidget::mouseDoubleClickEvent(QMouseEvent* e)
{
  mousePressEvent(e);
  if (m_activeTool == SelectTool && m_ltsManager->simulationActive() && m_ltsManager->selectedState())
  {
    m_ltsManager->simulateState(m_ltsManager->selectedState());
  }
}

void GLWidget::mouseMoveEvent(QMouseEvent* e)
{
  QPoint oldPosition = m_lastMousePosition;
  m_lastMousePosition = e->pos();

  if (!e->buttons() || !m_dragging)
  {
    e->ignore();
    return;
  }

  ArcballCamera* cam = dynamic_cast<ArcballCamera*>(m_scene.m_camera);

  if (m_activeTool == PanTool)
  {
    
    cam->applyTranslate(oldPosition, m_lastMousePosition, 1000);
    e->accept();
    repaint();
  }
  else if (m_activeTool == ZoomTool)
  {
    cam->applyZoom(oldPosition, m_lastMousePosition);
    e->accept();
    repaint();
  }
  else if (m_activeTool == RotateTool)
  {
    cam->applyDragRotate(oldPosition, m_lastMousePosition);
    mCRL2log(mcrl2::log::debug) << "Rotate..." << std::endl;
    e->accept();
    repaint();
  }
}

void GLWidget::onClusterStyleChanged()
{
  mCRL2log(mcrl2::log::debug) << "Cluster style changed" << std::endl;
  m_scene.rebuildScene();
}

void GLWidget::onDisplayOptionChanged()
{
  //m_scene.rebuildScene();
}
void GLWidget::onObjectsChanged()
{
  //m_scene.rebuildScene();
}

void GLWidget::rebuildScene()
{
  m_scene.rebuildScene();
}

void GLWidget::wheelEvent(QWheelEvent* e){
  ArcballCamera* cam = dynamic_cast<ArcballCamera*>(m_scene.m_camera);
  QPoint p1(0, 0);
  QPoint p2(0, -e->angleDelta().y());
  cam->applyZoom(p1, p2, 1.0005f);
  update();
};
