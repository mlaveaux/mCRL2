// Author(s): Ruben Vink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_LTSVIEW_GLWIDGET_H
#define MCRL2_LTSVIEW_GLWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLDebugLogger>
#include "scene.h"

/// \brief A Qt OpenGL widget which uses the GLScene interface to render the contents of a graph.
class GLWidget : public QOpenGLWidget
{
  Q_OBJECT
public:

  /**
   * @brief Constructor.
   * @param parent The parent widget for the user interface.
   */
  explicit GLWidget(Cluster* root, QWidget* parent = nullptr);

  /**
   * @brief Destructor.
   */
  ~GLWidget() override;

  /**
   * @brief Initialises the OpenGL context.
   */
  void initializeGL() override;

  /**
   * @brief Paints the OpenGL context.
   */
  void paintGL() override;

  /**
   * @brief Resize the OpenGL viewport.
   * @param width The new width for the viewport.
   * @param height The new height for the viewport.
   */
  void resizeGL(int width, int height) override;

  /**
   * @brief Processes mouse click events.
   * @param e The mouse event.
   */
  void mousePressEvent(QMouseEvent* e) override;

  /**
   * @brief Processes mouse click release events.
   * @param e The mouse event.
   */
  void mouseReleaseEvent(QMouseEvent* e) override;

  /**
   * @brief Processes mouse move events.
   * @param e The mouse event.
   */
  void mouseMoveEvent(QMouseEvent* e) override;

  /**
   * @brief Processes mouse scroll/wheel events.
   * @param e The mouse event.
   */
  void wheelEvent(QWheelEvent* e) override;

  /// TOOD: Concrete
  void update(Cluster* root);
private:
  GlLTSView::Scene m_scene;              ///< The GLScene which is used to render the contents.
  Cluster* m_root;                       ///< Root of cluster tree to be drawn on screen

  QOpenGLDebugLogger* m_ogl_logger;

  void logMessage(const QOpenGLDebugMessage& debugMessage);
};
#endif // MCRL2_LTSVIEW_GLWIDGET_H
