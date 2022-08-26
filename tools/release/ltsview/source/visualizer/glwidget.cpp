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

GLWidget::GLWidget(Cluster* root, QWidget* parent)
    : QOpenGLWidget(parent), m_scene(this, root)
{
  SceneGraph<GlLTSView::NodeData, GlLTSView::SceneData> sg;
  m_scene.m_camera = new ArcballCamera();
  m_scene.m_scenegraph = sg;
}

void GLWidget::update(Cluster* root)
{
  std::cout << "GLWidget::update(Cluster* root) called." << std::endl;
  makeCurrent();
  m_scene.m_clusterRoot = root;
  m_scene.rebuild();
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
    QPair<int, int> required(4, 3);
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
    {
      mCRL2log(mcrl2::log::verbose) << "Created an OpenGL " << version.first
                                    << "." << version.second << " context.\n";
    }
  }
  else
  {
    mcrl2::runtime_error("Context creation failed.");
  }

  // Enable real-time logging of OpenGL errors when the GL_KHR_debug extension
  // is available.
  m_ogl_logger = new QOpenGLDebugLogger(this);
  if (m_ogl_logger->initialize())
  {
    connect(m_ogl_logger, &QOpenGLDebugLogger::messageLogged, this,
            &GLWidget::logMessage);
    m_ogl_logger->startLogging();
  }
  else
  {
    mCRL2log(mcrl2::log::debug)
        << "QOpenGLDebugLogger initialisation failed.\n";
  }
  makeCurrent();
  m_scene.initialize();
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
  QPainter painter = QPainter(this);
  m_scene.render(painter);
};

void GLWidget::resizeGL(int width, int height){};

void GLWidget::mousePressEvent(QMouseEvent* e){};

void GLWidget::mouseReleaseEvent(QMouseEvent* e){};

void GLWidget::mouseMoveEvent(QMouseEvent* e){};

void GLWidget::wheelEvent(QWheelEvent* e){};
