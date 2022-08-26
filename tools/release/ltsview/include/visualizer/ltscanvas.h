// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef LTSCANVAS_H
#define LTSCANVAS_H

#include <QOpenGLFunctions_4_3_Core>
#include <QCursor>
#include <QOpenGLWidget>
#include <QImage>
#include <QMouseEvent>
#include <QPoint>
#include <QWheelEvent>
#include <QQuaternion>

#include <QLinkedList>
#include <QTime>

#include "visualizer.h"

class LtsCanvas : public QOpenGLWidget
{
  Q_OBJECT

  private:
    enum Tool { SelectTool, PanTool, ZoomTool, RotateTool };
    enum ObjectType { StateObject, ClusterObject };
    struct Selection { State* state; Cluster* cluster; };

  public:
    LtsCanvas(QWidget* parent, LtsManager* ltsManager, MarkManager* markManager);
    int viewWidth() const { return m_width; }
    int viewHeight() const { return m_height; }

  public slots:
    void resetView();
    void exportToText(QString filename) { m_visualizer->exportToText(filename.toStdString()); }

  protected slots:
    void clusterPositionsChanged();

  protected:
    void initializeGL() override;
    void resizeGL(int width, int height) override;
    void paintGL() override;
    void render(bool light);
    void determineActiveTool(QMouseEvent* event, bool useModifiers);
    void setActiveTool(Tool tool);
    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void mouseDoubleClickEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void wheelEvent(QWheelEvent* event);
    Selection selectObject(QPoint position);
    Selection parseSelection(GLuint* selectionBuffer, GLint items);
    void applyRotation(bool reverse = false);
    QVector3D getArcBallVector(int screenX, int screenY);

  public:
    QImage renderImage(int width, int height);
    void renderVectorStart() { emit renderingStarted(); }
    void renderVectorPage() { render(false); }
    void renderVectorFinish() { emit renderingFinished(); }

  signals:
    void renderingStarted();
    void renderingFinished();

  private:
    LtsManager* m_ltsManager;
    Visualizer* m_visualizer;

    QCursor m_selectCursor;
    QCursor m_panCursor;
    QCursor m_zoomCursor;
    QCursor m_rotateCursor;

    int m_width;
    int m_height;
    float m_baseDepth;
    float m_nearPlane;
    float m_farPlane;
    QVector3D m_position;
    QQuaternion m_rotation;
    Tool m_activeTool;
    bool m_dragging;
    QPoint m_lastMousePosition;
    QOpenGLFunctions_4_3_Core* m_ogl;
    QOpenGLContext* m_context;

    // graphics info variables
    QLinkedList<QTime> frame_times; 
    double measuring_time = 1; // [s]
  public:
    QString graphics_info;
};

#endif
