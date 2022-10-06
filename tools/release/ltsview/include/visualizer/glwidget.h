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
#include "testscene.h"
#include "scene.h"
#include "ltsmanager.h"

/// \brief A Qt OpenGL widget which uses the GLScene interface to render the contents of a graph.
// class GLWidget : public QOpenGLWidget, protected GlLTSView::Scene
class GLWidget : public QOpenGLWidget, protected Test::TScene
{
  Q_OBJECT  
  
private:
    enum Tool { SelectTool, PanTool, ZoomTool, RotateTool };
    enum ObjectType { StateObject, ClusterObject };
    struct Selection { State* state; Cluster* cluster; };

    void determineActiveTool(QMouseEvent* event, bool useModifiers);
    void setActiveTool(Tool tool);

    Selection selectObject(QPoint position);

    /// @brief Currently selected tool
    Tool m_activeTool;

    /// @brief different cursors for different operations
    QCursor m_selectCursor;
    QCursor m_panCursor;
    QCursor m_zoomCursor;
    QCursor m_rotateCursor;     
    
    bool m_dragging;            ///< Mouse currently held down and moving
    QPoint m_lastMousePosition; ///< Last mouse position

public:

  /**
   * @brief Constructor.
   * @param parent The parent widget for the user interface.
   */
  explicit GLWidget(Cluster* root, LtsManager* ltsManager, QWidget* parent = nullptr);

  /**
   * @brief Destructor.
   */
  ~GLWidget() override;

protected:
  /**
   * @brief Initialises the OpenGL context.
   */
  void initializeGL();

  /**
   * @brief Paints the OpenGL context.
   */
  void paintGL();

  /**
   * @brief Resize the OpenGL viewport.
   * @param width The new width for the viewport.
   * @param height The new height for the viewport.
   */
  void resizeGL(int width, int height);
  int m_width, m_height;
public:
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

  /**
   * @brief Processes double click mouse events. Calls single press event.
   * 
   * @param e The mouse event.
   */
  void mouseDoubleClickEvent(QMouseEvent* e);

  void setRoot(Cluster* root);
private:
  QOpenGLDebugLogger* m_ogl_logger;
  LtsManager* m_ltsManager;
  void logMessage(const QOpenGLDebugMessage& debugMessage);
};
#endif // MCRL2_LTSVIEW_GLWIDGET_H
