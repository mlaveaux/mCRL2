// Author(s): A.J. (Hannes) Pretorius
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file ./visualizer.cpp

#include <iostream> // only temporary for std::clog
#include <QOpenGLFramebufferObject>
#include <QOpenGLDebugLogger>

#include <QMessageBox>
#include "visualizer.h"
#include "mcrl2/utilities/logger.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

Visualizer::Visualizer(
  QWidget *parent,
  Graph *graph_)
  : QOpenGLWidget(parent),
    m_lastMouseEvent(std::make_unique<QMouseEvent>(QEvent::None, QPoint(0,0), Qt::NoButton, Qt::NoButton, Qt::NoModifier)),
    m_graph(graph_)
{
  setMinimumSize(10,10);
  setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  setFocusPolicy(Qt::ClickFocus);
  setUpdateBehavior(QOpenGLWidget::PartialUpdate);

  initMouse();

  texCharOK = false;
  texCushOK = false;

  geomChanged = true;
  dataChanged = true;

  showMenu = false;
}

void Visualizer::updateSelection() {
  const SelectionList selections = getSelection();
  if (selections.empty())
    handleSelection({});
  else
    handleSelection(selections.back());
}

Visualizer::SelectionList Visualizer::getSelection(qreal width, qreal height)
{
  makeCurrent();
  m_selectionBuffer->bind();

  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);

  constexpr std::size_t buffer_size = 4096;
  GLuint buffer[buffer_size];
  glSelectBuffer(buffer_size, buffer);
  glRenderMode(GL_SELECT);
  glInitNames();

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  const qreal x = m_lastMouseEvent->position().x() * devicePixelRatio();
  const qreal y = m_lastMouseEvent->position().y() * devicePixelRatio();
  gluPickMatrix(x, viewport[3] - y, width, height, viewport);
  GLdouble scale = std::min(widthC(), heightC());
  gluOrtho2D(-widthC() / scale, widthC() / scale, -heightC() / scale, heightC() / scale);

  glMatrixMode(GL_MODELVIEW);
  mark();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  GLint hits = glRenderMode(GL_RENDER);
  std::list<Selection> selections;
  GLuint* hit = buffer;
  while (hits--)
  {
    const size_t size = *hit++;
    hit += 2; // skip depth values
    selections.emplace_back(size);
    for (size_t i = 0; i < size; ++i)
    {
      selections.back()[i] = *hit++;
    }
  }

  m_selectionBuffer->release();
  doneCurrent();

  return selections;
}

void Visualizer::initializeGL() 
{
  if (mCRL2logEnabled(mcrl2::log::log_level_t::debug)) 
  {  
    QPair<int, int> version = format().version();
    qDebug() << "Created an OpenGL " << version.first
            << "." << version.second << " context.\n";
          
    // Enable real-time logging of OpenGL errors when the GL_KHR_debug extension
    // is available.
    m_logger = new QOpenGLDebugLogger(this);
    if (m_logger->initialize())
    {
      connect(m_logger, &QOpenGLDebugLogger::messageLogged, this,
              &Visualizer::logMessage);
      m_logger->startLogging();
    }
    else
    {
      qDebug() << "QOpenGLDebugLogger initialisation failed\n";
    }
  }
  m_selectionBuffer = std::make_unique<QOpenGLFramebufferObject>(1, 1);
  if (!m_selectionBuffer->isValid())
  {
    throw mcrl2::runtime_error("Failed to create framebuffer for selection handling.");
  }
}

void Visualizer::logMessage(const QOpenGLDebugMessage& debugMessage)
{
  qDebug() << "OpenGL: " << debugMessage.message() << "\n";
}

void Visualizer::paintGL()
{
  updateGeometry(); // TODO: move this
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  GLdouble scale = std::min(widthC(), heightC());
  gluOrtho2D(-widthC() / scale, widthC() / scale, -heightC() / scale, heightC() / scale);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glViewport(0, 0, widthC(), heightC());
  visualize();
}


QSizeF Visualizer::worldSize()
// Return viewport width and height in WORLD coordinates. The 
// viewport is set up such that the shortest side has
// length 2 in world coordinates.
{
  if (height() < width())
  {
    return QSizeF(width() / (double)height() * 2.0, 2.0);
  }
  else
  {
    return QSizeF(2.0, height() / (double)width() * 2.0);
  }
}


double Visualizer::pixelSize()
// Return distance in WORLD coordinates of 1 pixel.
{
  return worldSize().width() / (double)(width()*devicePixelRatio());
}


QPointF Visualizer::worldCoordinate(QPointF deviceCoordinate)
{
  QSizeF size = worldSize();
  double pixel = pixelSize();
  return QPointF(-0.5 * size.width() + deviceCoordinate.x() * pixel * devicePixelRatio(), 
                  0.5 * size.height() - deviceCoordinate.y() * pixel * devicePixelRatio());
}


// -- visualization functions ---------------------------------------




// -- event handlers ------------------------------------------------


void Visualizer::handleSizeEvent()
{
  geomChanged = true;
}

void Visualizer::handleMouseEvent(QMouseEvent* e)
{
  m_mouseDragReleased = false;
  if (!m_mouseDrag && e->buttons() != Qt::NoButton && e->type() == QEvent::MouseMove)
  {
    m_mouseDrag = true;
    m_mouseDragStart = m_lastMouseEvent->pos();
  }
  if (m_mouseDrag && e->buttons() == Qt::NoButton)
  {
    m_mouseDrag = false;
    m_mouseDragReleased = true;
  }

  m_lastMouseEvent = std::make_unique<QMouseEvent>(e->type(), e->pos(), e->globalPosition(), e->button(), e->buttons(), e->modifiers());
}

void Visualizer::handleKeyEvent(QKeyEvent* e)
{
  if (e->type() == QEvent::KeyPress)
  {
    m_lastKeyCode = Qt::Key(e->key());
  }
  else
  {
    m_lastKeyCode = Qt::Key_unknown;
  }
}


// -- protected utility functions -----------------------------------


void Visualizer::clear()
{
  VisUtils::clear(Qt::white);
}


void Visualizer::initMouse()
{
  m_lastMouseEvent = std::make_unique<QMouseEvent>(QEvent::None, QPoint(0,0), Qt::NoButton, Qt::NoButton, Qt::NoModifier);
  m_mouseDrag = false;
  m_mouseDragStart = QPoint(0,0);
}


void Visualizer::genCharTex()
{
  VisUtils::genCharTextures(
    texCharId,
    texChar);
  texCharOK = true;
}


void Visualizer::genCushTex()
{
  VisUtils::genCushTextures(
    texCushId,
    texCush);
  texCushOK = true;
}

void Visualizer::updateGeometry()
{
  m_canvas_width = width() * devicePixelRatio();
  m_canvas_height = height() * devicePixelRatio();
}

