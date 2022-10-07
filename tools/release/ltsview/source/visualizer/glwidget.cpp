// Author(s): Ruben Vink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <QtOpenGL>
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
    : QOpenGLWidget(parent), Test::TScene(root), m_ltsManager(ltsManager)
{
  SceneGraph<Test::NodeData, Test::SceneData> sg =
      SceneGraph<Test::NodeData, Test::SceneData>();
  m_camera = new ArcballCamera();
  m_scenegraph = sg;
}

void GLWidget::setRoot(Cluster* root)
{
  m_clusterRoot = root;
  rebuildScene();
}

GLWidget::~GLWidget(){};

/**
 * @brief Initialises the OpenGL context.
 */
void GLWidget::initializeGL()
{
  QOpenGLFunctions* f = QOpenGLContext::currentContext()->functions();
  f->glClearColor(0, 0, 0, 1);
  f->glEnable(GL_DEPTH_TEST);
  f->glEnable(GL_CULL_FACE);
  mCRL2log(mcrl2::log::debug) << "initializeGL called" << std::endl;
  initializeScene();
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
  renderScene();
  painter->endNativePainting();
};

void GLWidget::resizeGL(int width, int height)
{
  QOpenGLFunctions* f = QOpenGLContext::currentContext()->functions();
  f->glViewport(0, 0, width, height);
  m_width = width; m_height = height;
  m_camera->setViewport(width, height);
};

void GLWidget::mousePressEvent(QMouseEvent* e)
{
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

  ArcballCamera* cam = dynamic_cast<ArcballCamera*>(m_camera);

  if (m_activeTool == PanTool)
  {
    cam->applyTranslate(oldPosition, m_lastMousePosition);
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
    e->accept();
    repaint();
  }
}

void GLWidget::wheelEvent(QWheelEvent* e){
  ArcballCamera* cam = dynamic_cast<ArcballCamera*>(m_camera);
  QPoint p1(0, 0);
  QPoint p2(0, -e->angleDelta().y());
  cam->applyZoom(p1, p2, 1.0005f);
  update();
};
