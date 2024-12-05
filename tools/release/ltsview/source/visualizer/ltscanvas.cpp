// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "mcrl2/utilities/logger.h"
#include "visualizer/ltscanvas.h"
#include <tr/tr.h>

#include "mcrl2/gui/arcball.h"
#include "mcrl2/gui/glu.h"

#include "icons/zoom_cursor.xpm"
#include "icons/pan_cursor.xpm"
#include "icons/rotate_cursor.xpm"

// #include "scene.h"
#include "generic_visualizer/cameras/arcballcamera.h"
#include "state.h" 

#include "generic_visualizer/glutil/glutil.h"

#include <chrono>

LtsCanvas::LtsCanvas(QWidget* parent, LtsManager* ltsManager, MarkManager* markManager):
  QOpenGLWidget(parent),
  m_ltsManager(ltsManager),
  m_visualizer(new Visualizer(this, ltsManager, markManager)),
  m_width(0),
  m_height(0),
  m_dragging(false),
  frame_times(std::list<QTime>())
{
  m_selectCursor = QCursor(Qt::ArrowCursor);
  m_panCursor = QCursor(QPixmap(pan_cursor));
  m_zoomCursor = QCursor(QPixmap(zoom_cursor));
  m_rotateCursor = QCursor(QPixmap(rotate_cursor));
  

  connect(m_visualizer, SIGNAL(dirtied()), this, SLOT(update()));
  connect(m_ltsManager, SIGNAL(clusterPositionsChanged()), this, SLOT(clusterPositionsChanged()));
  connect(m_ltsManager, SIGNAL(ltsZoomed(LTS *)), this, SLOT(clusterPositionsChanged()));
  connect(m_ltsManager, SIGNAL(ltsZoomed(LTS *)), this, SLOT(resetView()));
  connect(&Settings::instance().backgroundColor, SIGNAL(changed(QColor)), this, SLOT(update()));
  connect(&Settings::instance().displayStates, SIGNAL(changed(bool)), this, SLOT(update()));
  connect(&Settings::instance().displayTransitions, SIGNAL(changed(bool)), this, SLOT(update()));
  connect(&Settings::instance().displayBackpointers, SIGNAL(changed(bool)), this, SLOT(update()));
  connect(&Settings::instance().displayWireframe, SIGNAL(changed(bool)), this, SLOT(update()));

  clusterPositionsChanged();
  resetView();
}

void LtsCanvas::resetView()
{
  m_position = QVector3D(0.0f, 0.0f, 0.0f);
  // structure will be drawn around the positive z-axis starting from the
  // origin, so we start with a rotation that makes the z-axis point downwards
  m_rotation = QQuaternion(1.0f, 1.0f, 0.0f, 0.0f).normalized();
  update();
}

void LtsCanvas::clusterPositionsChanged()
{
  float width;
  float height;
  m_visualizer->computeBoundsInfo(width, height);
  m_nearPlane = 1.0f;
  m_baseDepth = 0.866f * height + width + m_nearPlane;
  m_farPlane = 4 * m_baseDepth;
  update();
}

void LtsCanvas::determineActiveTool(QMouseEvent* event, bool useModifiers)
{
  if (useModifiers && (event->modifiers() & Qt::ControlModifier) && (event->modifiers() & Qt::ShiftModifier))
  {
    setActiveTool(ZoomTool);
  }
  else if (useModifiers && event->modifiers() & Qt::ControlModifier)
  {
    setActiveTool(PanTool);
  }
  else if ((event->buttons() & Qt::MiddleButton) || ((event->buttons() & Qt::LeftButton) && (event->buttons() & Qt::RightButton)))
  {
    setActiveTool(ZoomTool);
  }
  else if ((useModifiers && event->modifiers() & Qt::ShiftModifier) || event->buttons() & Qt::RightButton)
  {
    setActiveTool(RotateTool);
  }
  else
  {
    setActiveTool(SelectTool);
  }
}

void LtsCanvas::setActiveTool(Tool tool)
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

void LtsCanvas::initializeGL()
{
  mCRL2log(mcrl2::log::debug)
      << "Waarom zijn we in LTSCanvas::initializeGL?" << std::endl;
  return;
  GLfloat gray[] = { 0.35f, 0.35f, 0.35f, 1.0f };
  GLfloat light_pos[] = { 50.0f, 50.0f, 50.0f, 1.0f };
  glLightfv(GL_LIGHT0, GL_AMBIENT, gray);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, gray);
  glLightfv(GL_LIGHT0, GL_POSITION, light_pos);

  glEnable(GL_NORMALIZE);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glDepthMask(GL_TRUE);
  glBlendFunc(GL_ONE, GL_ZERO);
  glDisable(GL_BLEND);

  GLfloat light_col[] = { 0.2f, 0.2f, 0.2f };
  glMaterialfv(GL_FRONT, GL_SPECULAR, light_col);
  glMaterialf(GL_FRONT, GL_SHININESS, 8.0f);
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);

  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  glClearDepth(1.0);
}

void LtsCanvas::resizeGL(int width, int height)
{
  m_width = width;
  m_height = height;

  glViewport(0, 0, m_width, m_height);
}

void LtsCanvas::paintGL()
{
  //render(m_dragging);
  return;
  // GLint64 start;
  // GLint64 end;
  // glGetInteger64v(GL_TIMESTAMP, &start);
  auto start = std::chrono::high_resolution_clock::now();
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0f, (GLfloat)m_width / (GLfloat)m_height, m_nearPlane, m_farPlane);

  emit renderingStarted();
  render(m_dragging);
  emit renderingFinished();

  auto end = std::chrono::high_resolution_clock::now();
  // glGetInteger64v(GL_TIMESTAMP, &end);

  frame_times.emplace_back(QTime::currentTime());
  while(frame_times.back().msecsSinceStartOfDay() - frame_times.front().msecsSinceStartOfDay() > measuring_time * 1000) frame_times.pop_front();
  graphics_info = QString("Version: ").append((const char *)glGetString(GL_VERSION));
  graphics_info.append("\nAverage frame time: ").append(QString::number(measuring_time / frame_times.size()));
  graphics_info.append("\nFrame_times length: ").append(QString::number(frame_times.size()));
  graphics_info.append("\nLast frametime in msec: ").append(QString::number((end - start).count() / 1e6));
  
}

void LtsCanvas::render(bool light)
{
  return;
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);

  glClearColor(
    Settings::instance().backgroundColor.value().red() / 255.0,
    Settings::instance().backgroundColor.value().green() / 255.0,
    Settings::instance().backgroundColor.value().blue() / 255.0,
    1.0f
  );
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  if (!light || Settings::instance().navLighting.value())
  {
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
  }
  else
  {
    glDisable(GL_NORMALIZE);
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
  }

  if (!light || Settings::instance().navSmoothShading.value())
  {
    glShadeModel(GL_SMOOTH);
  }
  else
  {
    glShadeModel(GL_FLAT);
  }

  if (Settings::instance().displayWireframe.value())
  {
    glPolygonMode(GL_FRONT,GL_LINE);
  }
  else
  {
    glPolygonMode(GL_FRONT,GL_FILL);
  }

  // apply panning, zooming and rotating transformations
  glTranslatef(m_position.x(), m_position.y(), m_position.z() - m_baseDepth);

  mcrl2::gui::applyRotation(m_rotation);

  // translate along the z-axis to make the vertical center of the structure
  // end up in the current origin
  float halfHeight = m_visualizer->getHalfStructureHeight();
  glTranslatef(0.0f, 0.0f, -halfHeight);

  if (m_ltsManager->simulationActive())
  {
    glPushName(StateObject);
    m_visualizer->drawSimStates(
      m_ltsManager->simulationStateHistory(),
      m_ltsManager->currentSimulationState(),
      m_ltsManager->currentSimulationTransition()
    );
    glPopName();
  }
  if (!light || Settings::instance().navShowStates.value())
  {
    if (Settings::instance().displayStates.value())
    {
      // Identify that we are drawing states
      glPushName(StateObject);
      m_visualizer->drawStates(m_ltsManager->simulationActive());
      glPopName();
    }
  }

  // Disable lighting while drawing transitions, otherwise their colours
  // change with the viewpoint
  glDisable(GL_NORMALIZE);
  glDisable(GL_LIGHTING);
  glDisable(GL_LIGHT0);

  m_visualizer->drawTransitions(
    Settings::instance().displayTransitions.value() && (!light || Settings::instance().navShowTransitions.value()),
    Settings::instance().displayBackpointers.value() && (!light || Settings::instance().navShowBackpointers.value())
  );

  if (m_ltsManager->simulationActive())
  {
    // Draw transitions followed during simulation and the possible
    // transitions going out of the current state.
    // Identify that we are drawing selectable sim states in this mode.
    m_visualizer->drawSimTransitions(
      !light || Settings::instance().navShowTransitions.value(),
      !light || Settings::instance().navShowBackpointers.value(),
      m_ltsManager->simulationTransitionHistory(),
      m_ltsManager->simulationAvailableTransitions(),
      m_ltsManager->currentSimulationTransition()
    );
  }

  // Enable lighting again, if required
  if (!light || Settings::instance().navLighting.value())
  {
    glEnable(GL_NORMALIZE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
  }

  if (!light || Settings::instance().navTransparency.value())
  {
    // determine current viewpoint in world coordinates
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(0.0f, 0.0f, halfHeight);
    mcrl2::gui::applyRotation(m_rotation, /*reverse=*/true);
    glTranslatef(-m_position.x(), -m_position.y(), -m_position.z() + m_baseDepth);
    GLfloat matrix[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
    QVector3D viewpoint = QVector3D(matrix[12], matrix[13], matrix[14]);
    glPopMatrix();
    // sort clusters on distance to viewpoint
    m_visualizer->sortClusters(viewpoint);
  }

  // draw the structure
  glEnable(GL_BLEND);
  glDepthMask(GL_FALSE);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  // Identify that we are drawing clusters
  glPushName(ClusterObject);
  m_visualizer->drawStructure();
  glPopName();
  glDepthMask(GL_TRUE);
  glBlendFunc(GL_ONE, GL_ZERO);
  glDisable(GL_BLEND);

  glPopMatrix();

  GLenum error = glGetError();
  if (error != GL_NO_ERROR)
  {
    qWarning("OpenGL error: %d", error);
  }

  
}

void LtsCanvas::mousePressEvent(QMouseEvent* event)
{
  determineActiveTool(event, true);

  if (m_activeTool == SelectTool)
  {
    if (m_ltsManager->lts())
    {
      Selection selection = selectObject(event->pos());
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
    m_lastMousePosition = event->pos();
  }
  event->accept();
}

void LtsCanvas::mouseReleaseEvent(QMouseEvent* event)
{
  determineActiveTool(event, false);
  m_dragging = m_activeTool != SelectTool;
  event->accept();
  update();
}

void LtsCanvas::mouseDoubleClickEvent(QMouseEvent* event)
{
  mousePressEvent(event);
  if (m_activeTool == SelectTool && m_ltsManager->simulationActive() && m_ltsManager->selectedState())
  {
    m_ltsManager->simulateState(m_ltsManager->selectedState());
  }
}

void LtsCanvas::mouseMoveEvent(QMouseEvent* event)
{
  QPoint oldPosition = m_lastMousePosition;
  m_lastMousePosition = event->pos();

  if (!event->buttons() || !m_dragging)
  {
    event->ignore();
    return;
  }

  if (m_activeTool == PanTool)
  {
    m_position += QVector3D(
      -0.0015f * (m_baseDepth - m_position.z()) * (oldPosition.x() - event->position().x()),
       0.0015f * (m_baseDepth - m_position.z()) * (oldPosition.y() - event->position().y()),
       0.0f
    );
    event->accept();
    update();
  }
  else if (m_activeTool == ZoomTool)
  {
    m_position += QVector3D(
      0.0f,
      0.0f,
      0.01f * (m_baseDepth - m_position.z()) * (oldPosition.y() - event->position().y())
    );
    event->accept();
    update();
  }
  else if (m_activeTool == RotateTool)
  {
    m_rotation = mcrl2::gui::arcballRotation(oldPosition, m_lastMousePosition) * m_rotation;
    event->accept();
    update();
  }
}

void LtsCanvas::wheelEvent(QWheelEvent* event)
{
  m_position += QVector3D(0.0f, 0.0f, 0.001f * (m_baseDepth - m_position.z()) * event->angleDelta().y());
  event->accept();
  update();
}

LtsCanvas::Selection LtsCanvas::selectObject(QPoint position)
{
  return Selection();
  // In the worst case, all the objects in the frame are hit.
  // For each hit, the following needs to be recorded:
  // * The number of names on the stack
  // * The minimal depth of the hit object
  // * The maximal depth of the hit object
  // * The identifier of the type of object clicked
  // * Up to two numbers indicating the object selected

  makeCurrent();
  int maxItems = m_ltsManager->lts()->getNumClusters() + m_ltsManager->lts()->getNumStates();
  int bufferSize = maxItems * 6;
  GLuint* selectionBuffer = new GLuint[bufferSize];

  glSelectBuffer(bufferSize, selectionBuffer);
  glRenderMode(GL_SELECT);
  glInitNames();

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glViewport(0, 0, m_width, m_height);

  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  gluPickMatrix((GLdouble)position.x(), (GLdouble)(viewport[3] - position.y()), 3.0, 3.0, viewport);

  gluPerspective(60.0f, (GLfloat)m_width / (GLfloat)m_height, m_nearPlane, m_farPlane);

  render(false);

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glFlush();

  GLint hits = glRenderMode(GL_RENDER);
  Selection output = parseSelection(selectionBuffer, hits);
  delete[] selectionBuffer;
  return output;
}

LtsCanvas::Selection LtsCanvas::parseSelection(GLuint* selectionBuffer, GLint items)
{
  GLuint* buffer = selectionBuffer;
  // This method selects the object clicked.
  //
  // The buffer content per hit is encoded as follows:
  // buffer[0]: The number of names on the name stack at the moment of the hit
  // buffer[1]: The minimal depth of the object hit
  // buffer[2]: The maximal depth of the object. We are certainly not interested in this.
  // buffer[3]: The type of the object picked
  // buffer[4]: The first identifier of the object picked.
  // (buffer[5]: The second identifier of the object picked.)

  bool stateFound = false;
  float minimumStateDepth = -1;
  int stateID = -1;

  float clusterFound = false;
  float minimumClusterDepth = -1;
  int rank = -1;
  int position = -1;

  for (GLint i = 0; i < items; i++)
  {
    GLuint size = *buffer++;
    assert (size <= 3);

    float minimalDepth = ((float)(*buffer++)) / (float)0x7fffffff;
    buffer++;

    if (!size)
    {
      continue;
    }

    GLuint itemData[3] = {0, 0, 0};
    for (GLuint j = 0; j < size; j++)
    {
      itemData[j] = *buffer++;
    }

    if (itemData[0] == StateObject)
    {
      if (!stateFound || minimalDepth < minimumStateDepth)
      {
        stateFound = true;
        minimumStateDepth = minimalDepth;
        stateID = (int)itemData[1];
      }
    }
    else if (itemData[0] == ClusterObject)
    {
      if (!clusterFound || minimalDepth < minimumClusterDepth)
      {
        clusterFound = true;
        minimumClusterDepth = minimalDepth;
        rank = (int)itemData[1];
        position = (int)itemData[2];
      }
    }
  }

  Selection output;
  output.state = 0;
  output.cluster = 0;

  if (stateFound)
  {
    output.state = m_ltsManager->lts()->state(stateID);
  }
  else if (clusterFound)
  {
    output.cluster = m_ltsManager->lts()->cluster(rank, position);
  }

  return output;
}

QImage LtsCanvas::renderImage(int width, int height)
{
  mCRL2log(mcrl2::log::debug)
      << "Waraom zijn we in renderImage????" << std::endl;
  makeCurrent();
  glReadBuffer(GL_BACK);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);

  TRcontext* context = trNew();
  trTileSize(context, m_width, m_height, 0);
  trImageSize(context, width, height);

  unsigned char* imageData;
  try
  {
    imageData = new unsigned char[width * height * 4];
  }
  catch (const std::bad_alloc& exception)
  {
    trDelete(context);
    return QImage();
  }

  trImageBuffer(context, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8_REV, imageData);
  trRowOrder(context, TR_TOP_TO_BOTTOM);
  trPerspective(context, 60.0f, (GLfloat)m_width / (GLfloat)m_height, m_nearPlane, m_farPlane);

  emit renderingStarted();
  int more;
  do
  {
    trBeginTile(context);
    render(false);
    more = trEndTile(context);
  }
  while (more);
  trDelete(context);
  emit renderingFinished();

  // RESET VIEW
  glViewport(0, 0, m_width, m_height);

  QImage image = QImage(imageData, width, height, QImage::Format_ARGB32).mirrored();
  delete[] imageData;
  return image;
}
