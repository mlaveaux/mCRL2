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
    : QOpenGLWidget(parent), Test::TScene(root)
{
  SceneGraph<Test::NodeData, Test::SceneData> sg = SceneGraph<Test::NodeData, Test::SceneData>();
  m_camera = new ArcballCamera();
  m_scenegraph = sg;
}

void GLWidget::setRoot(Cluster* root){
  m_clusterRoot = root;
  rebuildScene();
}

GLWidget::~GLWidget(){};

/**
 * @brief Initialises the OpenGL context.
 */
void GLWidget::initializeGL()
{
  QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();
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
  QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();
  mCRL2log(mcrl2::log::debug) << "paintGL called." << std::endl;
  // Cull polygons that are facing away (back) from the camera, where their front is defined as counter clockwise by default, see glFrontFace, meaning that the
  // vertices that make up a triangle should be oriented counter clockwise to show the triangle.
  f->glDisable(GL_CULL_FACE);

  // Enable depth testing, so that we don't have to care too much about rendering in the right order.
  f->glEnable(GL_DEPTH_TEST);

  // Change the alpha blending function to make an alpha of 1 opaque and 0 fully transparent.
  f->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  f->glEnable(GL_BLEND);

  // Enable multisample antialiasing.
  f->glEnable(GL_MULTISAMPLE);

  f->glClearColor(Settings::instance().backgroundColor.value().redF(),
               Settings::instance().backgroundColor.value().greenF(),
               Settings::instance().backgroundColor.value().blueF(),
               1.0);
  f->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  renderScene();
};

void GLWidget::resizeGL(int width, int height){ 
  QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();
  f->glViewport(0, 0, width, height);
};

void GLWidget::mousePressEvent(QMouseEvent* e){};

void GLWidget::mouseReleaseEvent(QMouseEvent* e){};

void GLWidget::mouseMoveEvent(QMouseEvent* e){};

void GLWidget::wheelEvent(QWheelEvent* e){};
