// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "visualizer.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include "cluster.h"
#include "mathutils.h"
#include "state.h"
#include "transition.h"

#include <QtOpenGL>

#include "scene.h"
#include "arcballcamera.h"
#include "renderer.h"

#include "glcolor.h"

using namespace MathUtils;

#define SELECT_BLEND 0.3f
QColor SELECT_COLOR = QColor(255, 122, 0);

Visualizer::Visualizer(QObject *parent, LtsManager *ltsManager_, MarkManager* markManager_):
  QObject(parent),
  ltsManager(ltsManager_),
  markManager(markManager_),
  primitiveFactory()
{
  connect(&Settings::instance().stateSize, SIGNAL(changed(float)), this, SIGNAL(dirtied()));
  connect(&Settings::instance().clusterHeight, SIGNAL(changed(float)), this, SLOT(dirtyMatrices()));
  connect(&Settings::instance().branchRotation, SIGNAL(changed(int)), this, SLOT(dirtyMatrices()));
  connect(&Settings::instance().branchTilt, SIGNAL(changed(int)), this, SLOT(branchTiltChanged(int)));
  connect(&Settings::instance().quality, SIGNAL(changed(int)), this, SLOT(dirtyObjects()));
  connect(&Settings::instance().transparency, SIGNAL(changed(int)), this, SIGNAL(dirtied()));
  connect(&Settings::instance().stateColor, SIGNAL(changed(QColor)), this, SIGNAL(dirtied()));
  connect(&Settings::instance().downEdgeColor, SIGNAL(changed(QColor)), this, SIGNAL(dirtied()));
  connect(&Settings::instance().upEdgeColorFrom, SIGNAL(changed(QColor)), this, SIGNAL(dirtied()));
  connect(&Settings::instance().upEdgeColorTo, SIGNAL(changed(QColor)), this,
          SIGNAL(dirtied()));
  connect(&Settings::instance().longInterpolationUpEdge, SIGNAL(changed(bool)), this,
          SIGNAL(dirtied()));

  connect(&Settings::instance().markedColor, SIGNAL(changed(QColor)), this, SLOT(dirtyColorsMark()));
  connect(&Settings::instance().clusterColorTop, SIGNAL(changed(QColor)), this, SLOT(dirtyColorsNoMark()));
  connect(&Settings::instance().clusterColorBottom, SIGNAL(changed(QColor)), this, SLOT(dirtyColorsNoMark()));
  connect(&Settings::instance().longInterpolationCluster, SIGNAL(changed(bool)), this, SLOT(dirtyColorsNoMark()));
  connect(&Settings::instance().simPrevColor, SIGNAL(changed(QColor)), this, SIGNAL(dirtied()));
  connect(&Settings::instance().simCurrColor, SIGNAL(changed(QColor)), this, SIGNAL(dirtied()));
  connect(&Settings::instance().simSelColor, SIGNAL(changed(QColor)), this, SIGNAL(dirtied()));
  connect(&Settings::instance().simPosColor, SIGNAL(changed(QColor)), this, SIGNAL(dirtied()));
  connect(&Settings::instance().clusterVisStyleTubes, SIGNAL(changed(bool)), this, SLOT(dirtyObjects()));

  connect(ltsManager, SIGNAL(ltsChanged(LTS *)), this, SLOT(setClusterHeight()));
  connect(ltsManager, SIGNAL(clustersChanged()), this, SLOT(dirtyObjects()));
  connect(ltsManager, SIGNAL(ltsZoomed(LTS *)), this, SLOT(dirtyObjects()));
  connect(ltsManager, SIGNAL(clusterPositionsChanged()), this, SLOT(dirtyObjects()));
  connect(ltsManager, SIGNAL(statePositionsChanged()), this, SLOT(dirtyPositions()));
  connect(ltsManager, SIGNAL(selectionChanged()), this, SLOT(dirtyColors()));
  connect(ltsManager, SIGNAL(simulationChanged()), this, SLOT(dirtyColors()));

  connect(markManager, SIGNAL(marksChanged()), this, SLOT(dirtyColors()));

  sin_obt = float(sin(deg_to_rad(Settings::instance().branchTilt.value())));
  cos_obt = float(cos(deg_to_rad(Settings::instance().branchTilt.value())));

  update_objects = true;
  update_matrices = false;
  update_colors = false;
  update_positions = false;
}

float Visualizer::getHalfStructureHeight() const
{
  if (!ltsManager->lts())
  {
    return 0.0f;
  }
  return Settings::instance().clusterHeight.value() * (ltsManager->lts()->getNumRanks() - 1) / 2.0f;
}

void Visualizer::setClusterHeight()
{
  float ratio = ltsManager->lts()->getInitialState()->getCluster()->getBCRadius() /
                ltsManager->lts()->getInitialState()->getCluster()->getBCHeight();
  Settings::instance().clusterHeight.setValue(std::max(4,round_to_int(40.0f * ratio)) / 10.0f);
  dirtyObjects();
}

void Visualizer::branchTiltChanged(int value)
{
  sin_obt = float(sin(deg_to_rad(value)));
  cos_obt = float(cos(deg_to_rad(value)));
  dirtyMatrices();
}

void Visualizer::computeBoundsInfo(float& bcw,float& bch)
{
  bcw = 0.0f;
  bch = 0.0f;
  if (ltsManager->lts())
  {
    computeSubtreeBounds(ltsManager->lts()->getInitialState()->getCluster(), bcw, bch);
  }
}

void Visualizer::computeSubtreeBounds(Cluster* root, float& bw, float& bh)
{
  // compute the bounding cylinder of the structure.
  if (!root->hasDescendants())
  {
    bw = root->getTopRadius();
    bh = 2.0f * root->getTopRadius();
  }
  else
  {
    Cluster* desc;
    int i;
    for (i = 0; i < root->getNumDescendants(); ++i)
    {
      desc = root->getDescendant(i);
      if (desc != NULL)
      {
        if (desc->isCentered())
        {
          // descendant is centered
          float dw = 0.0f;
          float dh = 0.0f;
          computeSubtreeBounds(desc,dw,dh);
          bw = std::max(bw,dw);
          bh = std::max(bh,dh);
        }
        else
        {
          float dw = 0.0f;
          float dh = 0.0f;
          computeSubtreeBounds(desc,dw,dh);
          bw = std::max(bw,root->getBaseRadius() + dh*sin_obt + dw*cos_obt);
          bh = std::max(bh,dh*cos_obt + dw*sin_obt);
        }
      }
    }
    bw = std::max(bw,root->getTopRadius());
    bh += Settings::instance().clusterHeight.value();
  }
}

// ------------- STRUCTURE -----------------------------------------------------

void Visualizer::drawStructure()
{
  if (!ltsManager->lts())
  {
    return;
  }
  if (update_objects || update_matrices)
  {
    traverseTree();
  }
  if (update_colors)
  {
    updateColors();
  }
  visObjectFactory.drawObjects(&primitiveFactory,(int)((100 - Settings::instance().transparency.value()) * 2.55f),
                                markManager->stateMatchStyle() == MATCH_MULTI);
}

void Visualizer::traverseTree()
{
  if (update_objects)
  {
    visObjectFactory.clear();
    update_colors = true;
    QOpenGLWidget oglw;
    SceneGraph<GlLTSView::NodeData, GlLTSView::SceneData> sg;
    ArcballCamera camera;
    LTSRenderer renderer;
    GlLTSView::Scene scene = GlLTSView::Scene(oglw, sg, camera, ltsManager->lts()->getInitialState()->getCluster());
    scene.rebuild();
  }
  glPushMatrix();
  glLoadIdentity();
  Cluster *root = ltsManager->lts()->getInitialState()->getCluster();
  if (Settings::instance().clusterVisStyleTubes.value())
  {
    traverseTreeT(root, true, 0);
  }
  else
  {
    traverseTreeC(root, true, 0);
  }
  glPopMatrix();
  update_objects = false;
  update_matrices = false;
  update_positions = true;
}

void Visualizer::traverseTreeC(Cluster* root,bool topClosed,int rot)
{
  if (!root->hasDescendants())
  {
    float r = root->getTopRadius();
    glPushMatrix();
    glScalef(r,r,r);
    std::vector<int> ids;
    ids.push_back(root->getRank());
    ids.push_back(root->getPositionInRank());
    if (update_objects)
    {
      root->setVisObject(visObjectFactory.makeObject(
                           primitiveFactory.makeSphere(),ids));
    }
    else
    {
      visObjectFactory.updateObjectMatrix(root->getVisObject());
    }
    glPopMatrix();
  }
  else
  {
    int drot = rot + Settings::instance().branchRotation.value();
    if (drot >= 360)
    {
      drot -= 360;
    }
    if (update_objects)
    {
      root->clearBranchVisObjects();
    }

    glTranslatef(0.0f,0.0f,Settings::instance().clusterHeight.value());
    for (int i = 0; i < root->getNumDescendants(); ++i)
    {
      Cluster* desc = root->getDescendant(i);
      if (desc != NULL)
      {
        if (desc->isCentered())
        {
          if (root->getNumDescendants() > 1)
          {
            traverseTreeC(desc,false,drot);
          }
          else
          {
            traverseTreeC(desc,false,rot);
          }
        }
        else
        {
          glRotatef(-desc->getPosition() - rot, 0.0f, 0.0f, 1.0f);
          glTranslatef(root->getBaseRadius(), 0.0f, 0.0f);
          glRotatef(Settings::instance().branchTilt.value(), 0.0f, 1.0f, 0.0f);
          traverseTreeC(desc, true, drot);
          glRotatef(-Settings::instance().branchTilt.value(), 0.0f, 1.0f, 0.0f);
          glTranslatef(-root->getBaseRadius(), 0.0f, 0.0f);
          glRotatef(desc->getPosition() + rot, 0.0f, 0.0f, 1.0f);
        }
      }
    }
    glTranslatef(0.0f,0.0f,-Settings::instance().clusterHeight.value());

    float r = root->getBaseRadius() / root->getTopRadius();
    glPushMatrix();
    glTranslatef(0.0f,0.0f,0.5f*Settings::instance().clusterHeight.value());
    if (r > 1.0f)
    {
      r = 1.0f / r;
      glRotatef(180.0f,1.0f,0.0f,0.0f);
      glScalef(root->getBaseRadius(),root->getBaseRadius(),
               Settings::instance().clusterHeight.value());

      std::vector<int> ids;
      ids.push_back(root->getRank());
      ids.push_back(root->getPositionInRank());

      if (update_objects)
      {
        root->setVisObject(visObjectFactory.makeObject(
                             primitiveFactory.makeTruncatedCone(r,topClosed,
                                 root->getNumDescendants() > 1 || root->hasSeveredDescendants()),
                             ids));
      }
      else
      {
        visObjectFactory.updateObjectMatrix(root->getVisObject());
      }
      glPopName();
      glPopName();
    }
    else
    {
      glScalef(root->getTopRadius(),root->getTopRadius(),Settings::instance().clusterHeight.value());

      std::vector<int> ids;
      ids.push_back(root->getRank());
      ids.push_back(root->getPositionInRank());

      if (update_objects)
      {
        root->setVisObject(visObjectFactory.makeObject(
                             primitiveFactory.makeTruncatedCone(r,
                                 root->getNumDescendants() > 1 || root->hasSeveredDescendants(),
                                 topClosed), ids));
      }
      else
      {
        visObjectFactory.updateObjectMatrix(root->getVisObject());
      }
      glPopName();
      glPopName();
    }
    glPopMatrix();
  }
}

void Visualizer::traverseTreeT(Cluster* root, bool topClosed, int rot)
{
  if (!root->hasDescendants())
  {
    // root has no descendants; so draw it as a hemispheroid
    glPushMatrix();
    glScalef(root->getTopRadius(),root->getTopRadius(),
             std::min(root->getTopRadius(),Settings::instance().clusterHeight.value()));
    if (update_objects)
    {
      std::vector<int> ids;
      ids.push_back(root->getRank());
      ids.push_back(root->getPositionInRank());
      if (root == ltsManager->lts()->getInitialState()->getCluster())
      {
        // exception: draw root as a sphere if it is the initial cluster
        root->setVisObject(visObjectFactory.makeObject(
                             primitiveFactory.makeSphere(), ids));
      }
      else
      {
        root->setVisObject(visObjectFactory.makeObject(
                             primitiveFactory.makeHemisphere(),ids));
      }
    }
    else
    {
      visObjectFactory.updateObjectMatrix(root->getVisObject());
    }
    glPopMatrix();
  }
  else
  {
    int drot = rot + Settings::instance().branchRotation.value();
    if (drot >= 360)
    {
      drot -= 360;
    }
    float baserad = 0.0f;
    if (update_objects)
    {
      root->clearBranchVisObjects();
    }

    for (int i = 0; i < root->getNumDescendants(); ++i)
    {
      Cluster* desc = root->getDescendant(i);
      if (desc != NULL)
      {
        if (desc->isCentered())
        {
          if (update_objects)
          {
            root->addBranchVisObject(-1);
          }
          baserad = desc->getTopRadius();
          glTranslatef(0.0f,0.0f,Settings::instance().clusterHeight.value());
          if (root->getNumDescendants() > 1)
          {
            traverseTreeT(desc,false,drot);
          }
          else
          {
            traverseTreeT(desc,false,rot);
          }
          glTranslatef(0.0f,0.0f,-Settings::instance().clusterHeight.value());
        }
        else
        {
          // make the connecting cone
          float d_rad = root->getBaseRadius() - root->getTopRadius();
          float sz = sqrt(Settings::instance().clusterHeight.value() *
                          Settings::instance().clusterHeight.value() + d_rad * d_rad);
          float alpha, sign;
          if (d_rad < 0.0f)
          {
            sign = -1.0f;
            alpha = atan(Settings::instance().clusterHeight.value() / -d_rad);
          }
          else
          {
            sign = 1.0f;
            if (d_rad > 0.0f)
            {
              alpha = atan(Settings::instance().clusterHeight.value() / d_rad);
            }
            else
            {
              alpha = 0.5f * PI;
            }
          }
          glRotatef(-desc->getPosition()-rot,0.0f,0.0f,1.0f);
          glPushMatrix();
          glTranslatef(root->getTopRadius(),0.0f,0.0f);
          glRotatef(sign*(90.0f-rad_to_deg(alpha)),0.0f,1.0f,0.0f);
          glScalef(sz,sz,sz);
          if (update_objects)
          {
            std::vector<int> ids;
            ids.push_back(root->getRank());
            ids.push_back(root->getPositionInRank());
            root->addBranchVisObject(visObjectFactory.makeObject(
                                       primitiveFactory.makeObliqueCone(alpha,
                                           desc->getTopRadius()/sz,sign),ids));
          }
          else
          {
            visObjectFactory.updateObjectMatrix(root->getBranchVisObject(i));
          }
          glPopMatrix();

          // recurse into the subtree
          glTranslatef(root->getBaseRadius(),0.0f,
                       Settings::instance().clusterHeight.value());
          glRotatef(Settings::instance().branchTilt.value(),0.0f,1.0f,0.0f);
          traverseTreeT(desc,false,drot);
          glRotatef(-Settings::instance().branchTilt.value(),0.0f,1.0f,0.0f);
          glTranslatef(-root->getBaseRadius(),0.0f,
                       -Settings::instance().clusterHeight.value());
          glRotatef(desc->getPosition()+rot,0.0f,0.0f,1.0f);
        }
      }
    }

    if (baserad <= 0.0f && !root->hasSeveredDescendants())
    {
      // root has no centered descendant, so draw it as a hemispheroid
      glPushMatrix();
      glScalef(root->getTopRadius(),root->getTopRadius(),
               std::min(root->getTopRadius(),Settings::instance().clusterHeight.value()));
      if (update_objects)
      {
        std::vector<int> ids;
        ids.push_back(root->getRank());
        ids.push_back(root->getPositionInRank());
        if (root == ltsManager->lts()->getInitialState()->getCluster())
        {
          // exception: draw root as a sphere if it is the initial cluster
          root->setVisObject(visObjectFactory.makeObject(
                               primitiveFactory.makeSphere(), ids));
        }
        else
        {
          root->setVisObject(visObjectFactory.makeObject(
                               primitiveFactory.makeHemisphere(), ids));
        }
      }
      else
      {
        visObjectFactory.updateObjectMatrix(root->getVisObject());
      }
      glPopMatrix();
    }
    else
    {
      // root has centered descendant; so draw it as a truncated cone
      float r = baserad / root->getTopRadius();
      if (root->hasSeveredDescendants())
      {
        r = 1.0f;
      }
      glPushMatrix();
      glTranslatef(0.0f,0.0f,0.5f*Settings::instance().clusterHeight.value());
      if (r > 1.0f)
      {
        r = 1.0f / r;
        glRotatef(180.0f,1.0f,0.0f,0.0f);
        glScalef(baserad,baserad,Settings::instance().clusterHeight.value());
        if (update_objects)
        {
          std::vector<int> ids;
          ids.push_back(root->getRank());
          ids.push_back(root->getPositionInRank());
          root->setVisObject(visObjectFactory.makeObject(
                               primitiveFactory.makeTruncatedCone(r,topClosed,
                                   root->hasSeveredDescendants()), ids));
        }
        else
        {
          visObjectFactory.updateObjectMatrix(root->getVisObject());
        }
      }
      else
      {
        glScalef(root->getTopRadius(),root->getTopRadius(),
                 Settings::instance().clusterHeight.value());
        if (update_objects)
        {
          std::vector<int> ids;
          ids.push_back(root->getRank());
          ids.push_back(root->getPositionInRank());
          root->setVisObject(visObjectFactory.makeObject(
                               primitiveFactory.makeTruncatedCone(r,
                                   root->hasSeveredDescendants(),topClosed),ids));
        }
        else
        {
          visObjectFactory.updateObjectMatrix(root->getVisObject());
        }
      }
      glPopMatrix();
    }
  }
}

void Visualizer::updateColors()
{
  Cluster* cl;
  QColor c;
  if (markManager->markStyle() == NO_MARKS)
  {
    QColor from = Settings::instance().clusterColorTop.value();
    QColor to = Settings::instance().clusterColorBottom.value();

    float hueFrom = from.hueF();
    float hueTo = to.hueF();
    bool useLongInterpolationCluster = Settings::instance().longInterpolationCluster.value();

    float saturationFrom = from.saturationF();
    float saturationDelta = to.saturationF() - saturationFrom;
    float valueFrom = from.valueF();
    float valueDelta = to.valueF() - valueFrom;

    int ranks = ltsManager->lts()->getMaxRanks();

    for (Cluster_iterator ci = ltsManager->lts()->getClusterIterator(); !ci.is_end();
         ++ci)
    {
      cl = *ci;

      if (ranks == 1)
      {
        c = from;
      }
      else
      {
        float t = cl->getRank() / (float)(ranks);
        if (useLongInterpolationCluster){
          c = GlUtil::Color::lerp<GlUtil::Color::InterpolateMode::LONG, GlUtil::Color::ColorMode::HSV>(from, to, t);
        }else{
          c = GlUtil::Color::lerp<GlUtil::Color::InterpolateMode::SHORT,
                                  GlUtil::Color::ColorMode::HSV>(from, to, t);
        }
      }

      if (cl == ltsManager->selectedCluster())
      {
        c = GlUtil::Color::blend<GlUtil::Color::ColorMode::RGB>(c, SELECT_COLOR, SELECT_BLEND, 0);
      }
      // set color of cluster cl
      visObjectFactory.updateObjectColor(cl->getVisObject(),c);
      // and its branches
      for (int i = 0; i < cl->getNumBranchVisObjects(); ++i)
      {
        if (cl->getBranchVisObject(i) != -1)
        {
          visObjectFactory.updateObjectColor(
            cl->getBranchVisObject(i),c);
        }
      }
    }
  }
  else // markManager->markStyle() != NO_MARKS
  {
    for (Cluster_iterator ci = ltsManager->lts()->getClusterIterator(); !ci.is_end();
         ++ci)
    {
      cl = *ci;
      std::vector<QColor> rule_colours;

      c = QColor(255, 255, 255);
      if (markManager->isMarked(cl))
      {
        if (markManager->markStyle() == MARK_STATES &&
            markManager->stateMatchStyle() == MATCH_MULTI)
        {
          QList<QColor> colors = markManager->markColors(cl);
          rule_colours.reserve(colors.size());
          rule_colours.insert(rule_colours.end(), colors.begin(), colors.end());
        }
        else
        {
          c = Settings::instance().markedColor.value();
        }
      }

      if (cl == ltsManager->selectedCluster())
      {
        c = GlUtil::Color::blend<GlUtil::Color::ColorMode::RGB>(c, SELECT_COLOR, SELECT_BLEND, 0);
      }

      visObjectFactory.updateObjectColor(cl->getVisObject(),c);
      visObjectFactory.updateObjectTexture(cl->getVisObject(),
                                            rule_colours);

      for (int i = 0; i < cl->getNumBranchVisObjects(); ++i)
      {
        if (cl->getBranchVisObject(i) != -1)
        {
          visObjectFactory.updateObjectColor(
            cl->getBranchVisObject(i),c);
          visObjectFactory.updateObjectTexture(
            cl->getBranchVisObject(i),rule_colours);
        }
      }
    }
  }
  update_colors = false;
}

void Visualizer::sortClusters(const QVector3D& viewpoint)
{
  visObjectFactory.sortObjects(viewpoint);
}

// ------------- STATES --------------------------------------------------------

void Visualizer::drawStates(bool simulating)
{
  if (!ltsManager->lts())
  {
    return;
  }

  // Compute absolute positions of nodes, if necessary
  computeAbsPos();
  drawStates(ltsManager->lts()->getInitialState()->getCluster(), simulating);
}

void Visualizer::drawSimStates(QList<State*> historicStates,
                               State* currState, Transition* chosenTrans)
{
  if (!ltsManager->lts())
  {
    return;
  }

  // Compute absolute positions of nodes, if necessary
  computeAbsPos();

  float ns = Settings::instance().stateSize.value();
  QColor c;
  QVector3D p;

  std::set<State*> drawnStates;

  // Draw the current state.
  //QColor hisStateColor = settings->simStateColor.value();
  if (ltsManager->lts()->getZoomLevel() == currState->getZoomLevel())
  {
    if (markManager->isMarked(currState))
    {
      c = Settings::instance().markedColor.value();
    }
    else
    {
      c = Settings::instance().simCurrColor.value();
    }

    if (currState == ltsManager->selectedState())
    {
        c = GlUtil::Color::blend<GlUtil::Color::ColorMode::RGB>(c, SELECT_COLOR, SELECT_BLEND, 0);
    }

    glColor4ub(c.red(), c.green(), c.blue(), 255);
    QVector3D p = currState->getPositionAbs();

    glPushName(currState->getID());

    glPushMatrix();
    glTranslatef(p.x(), p.y(), p.z());

    // Make the current state a bit larger, to make it easier to find it in the
    // simulation
    glScalef(1.5 *ns, 1.5* ns, 1.5 * ns);
    primitiveFactory.drawSimpleSphere();
    glPopMatrix();

    glPopName();
    drawnStates.insert(currState);
  }

  // Draw the future states
  State* endState;
  for (int i = 0; i < currState->getNumOutTransitions(); ++i)
  {
    endState = currState->getOutTransition(i)->getEndState();

    if (ltsManager->lts()->getZoomLevel() == endState->getZoomLevel()
        && drawnStates.find(endState) == drawnStates.end())
    {
      glPushName(endState->getID());

      if (markManager->isMarked(endState))
      {
        c = Settings::instance().markedColor.value();
      }
      else if (currState->getOutTransition(i) == chosenTrans)
      {
        c = Settings::instance().simSelColor.value();
      }
      else
      {
        c = Settings::instance().simPosColor.value();
      }

      if (endState == ltsManager->selectedState())
      {
        c = GlUtil::Color::blend<GlUtil::Color::ColorMode::RGB>(c, SELECT_COLOR, SELECT_BLEND, 0);
      }

      glColor4ub(c.red(), c.green(), c.blue(), 255);
      p = endState->getPositionAbs();
      glPushMatrix();
      glTranslatef(p.x(), p.y(), p.z());
      if (currState->getOutTransition(i) == chosenTrans)
      {
        glScalef(1.5 * ns, 1.5 * ns, 1.5* ns);
      }
      else
      {
        glScalef(ns, ns, ns);
      }
      primitiveFactory.drawSimpleSphere();
      glPopMatrix();

      glPopName();
      drawnStates.insert(endState);
    }
  }

  // Draw previous states of the simulation, in the colour specified in the
  // settings (default: white)
  State* s;
  for (int i = 0; i < historicStates.size(); ++i)
  {
    s = historicStates[i];

    if (ltsManager->lts()->getZoomLevel() == s->getZoomLevel()
        && drawnStates.find(s) == drawnStates.end())
    {
      if (markManager->isMarked(s))
      {
        c = Settings::instance().markedColor.value();
      }
      else
      {
        c = Settings::instance().simPrevColor.value();
      }

      if (s == ltsManager->selectedState())
      {
        c = GlUtil::Color::blend<GlUtil::Color::ColorMode::RGB>(c, SELECT_COLOR, SELECT_BLEND, 0);
      }

      glColor4ub(c.red(), c.green(), c.blue(), 255);

      glPushName(s->getID());

      p = s->getPositionAbs();
      glPushMatrix();
      glTranslatef(p.x(), p.y(), p.z());
      glScalef(ns, ns, ns);
      primitiveFactory.drawSimpleSphere();
      glPopMatrix();

      glPopName();
      drawnStates.insert(s);
    }
  }
}

void Visualizer::computeAbsPos()
{
  if (update_objects || update_matrices)
  {
    traverseTree();
  }
  if (update_positions)
  {
    glPushMatrix();
    glLoadIdentity();
    computeStateAbsPos(ltsManager->lts()->getInitialState()->getCluster(), 0);
    glPopMatrix();
    update_positions = false;
  }
}

void Visualizer::computeStateAbsPos(Cluster* root, int rot)
{
// Does a DFS on the clusters to calculate the `absolute' positions of their
// states, taking the position of the initial state as (0,0,0).
// Secondly, it also assigns the incoming and outgoing control points for the
// back pointers to each state.

  if (root->getState(0)->getZoomLevel() == ltsManager->lts()->getZoomLevel())
  {
    // compute position of each state in this cluster
    for (int i = 0; i < root->getNumStates(); ++i)
    {
      State* s = root->getState(i);
      if (s->isCentered())
      {
        float M[16];

        // store the position of this state
        glGetFloatv(GL_MODELVIEW_MATRIX, M);
        QVector3D p1 = QVector3D(M[12], M[13], M[14]);
        s->setPositionAbs(p1);

        // The outgoing vector of the state lies settings->clusterHeight.value() above the state.
        glPushMatrix();
        glTranslatef(0.0f, 0.0f, -2 * Settings::instance().clusterHeight.value());
        glGetFloatv(GL_MODELVIEW_MATRIX, M);
        QVector3D p2 = QVector3D(M[12], M[13], M[14]);
        s->setIncomingControl(p2);

        // The incoming vector of the state lies settings->clusterHeight.value() beneath the state.
        glTranslatef(0.0f, 0.0f, 4 * Settings::instance().clusterHeight.value());
        glGetFloatv(GL_MODELVIEW_MATRIX, M);
        QVector3D p3 = QVector3D(M[12], M[13], M[14]);
        s->setOutgoingControl(p3);

        glPopMatrix();
      }
      else
      {
        float M[16];
        glPushMatrix();
        glRotatef(-s->getPositionAngle(), 0.0f, 0.0f, 1.0f);

        glTranslatef(s->getPositionRadius(), 0.0f, 0.0f);
        glGetFloatv(GL_MODELVIEW_MATRIX, M);
        QVector3D p1 = QVector3D(M[12], M[13], M[14]);
        s->setPositionAbs(p1);
        glTranslatef(-s->getPositionRadius(), 0.0f, 0.0f);

        // The outgoing vector of the state points out of the cluster, in the
        // direction the state itself is positioned. Furthermore, it points
        // settings->clusterHeight.value() up.
        glTranslatef(root->getTopRadius() * 3, 0.0f, -Settings::instance().clusterHeight.value());
        glGetFloatv(GL_MODELVIEW_MATRIX, M);
        QVector3D p2 = QVector3D(M[12], M[13], M[14]);
        s->setIncomingControl(p2);

        glTranslatef(0.0f, 0.0f, 2 * Settings::instance().clusterHeight.value());
        glGetFloatv(GL_MODELVIEW_MATRIX, M);
        QVector3D p3 = QVector3D(M[12], M[13], M[14]);
        s->setOutgoingControl(p3);
        glPopMatrix();
      }

      float M[16];
      glPushMatrix();
      glRotatef(-s->getPositionAngle(), 0.0f, 0.0f, 1.0f);
      glTranslatef(s->getPositionRadius(), 0.0f, 0.0f);

      glTranslatef(Settings::instance().stateSize.value() * 5.0f,
                   Settings::instance().stateSize.value() * 5.0f,
                   0.0f);
      glGetFloatv(GL_MODELVIEW_MATRIX, M);
      QVector3D p = QVector3D(M[12], M[13], M[14]);
      s->setLoopControl1(p);

      glTranslatef(0.0f,
                   -Settings::instance().stateSize.value() * 10.0f,
                   0.0f);
      glGetFloatv(GL_MODELVIEW_MATRIX, M);
      p = QVector3D(M[12], M[13], M[14]);
      s->setLoopControl2(p);

      glPopMatrix();
    }

    // recurse into the descendants
    int drot = rot + Settings::instance().branchRotation.value();
    if (drot >= 360)
    {
      drot -= 360;
    }
    glTranslatef(0.0f,0.0f,Settings::instance().clusterHeight.value());

    for (int i = 0; i < root->getNumDescendants(); ++i)
    {
      Cluster* desc = root->getDescendant(i);
      if (desc != NULL)
      {
        if (desc->isCentered())
        {
          computeStateAbsPos(desc,(root->getNumDescendants()>1)?drot:rot);
        }
        else
        {
          glRotatef(-desc->getPosition()-rot,0.0f,0.0f,1.0f);
          glTranslatef(root->getBaseRadius(),0.0f,0.0f);
          glRotatef(Settings::instance().branchTilt.value(),0.0f,1.0f,0.0f);
          computeStateAbsPos(desc,drot);
          glRotatef(-Settings::instance().branchTilt.value(),0.0f,1.0f,0.0f);
          glTranslatef(-root->getBaseRadius(),0.0f,0.0f);
          glRotatef(desc->getPosition()+rot,0.0f,0.0f,1.0f);
        }
      }
    }

    glTranslatef(0.0f,0.0f,-Settings::instance().clusterHeight.value());
  }
}

void Visualizer::drawStates(Cluster* root, bool simulating)
{
  if (ltsManager->lts()->getZoomLevel() == root->getState(0)->getZoomLevel())
  {
    float ns = Settings::instance().stateSize.value();
    for (int i = 0; i < root->getNumStates(); ++i)
    {
      State* s = root->getState(i);

      if (!(simulating && s->isSimulated()))
      {
        QColor c;
        GLuint texName;

        if (!markManager->isMarked(s))
        {
          c = Settings::instance().stateColor.value();
        }
        else
        {
          if (markManager->stateMatchStyle() != MATCH_MULTI)
          {
            c = Settings::instance().markedColor.value();
          }
          else
          {
            c = Settings::instance().stateColor.value();

            QList<QColor> colors = markManager->markColors(s);

            int n_colours = 1;
            while (n_colours < colors.size())
            {
              n_colours = n_colours << 1;
            }
            GLubyte* texture = (GLubyte*)malloc(4 * n_colours *
                                                sizeof(GLubyte));

            for (int i = 0; i < n_colours; ++i)
            {
              QColor c1 = colors[i % colors.size()];
              texture[4*i]   = (GLubyte) c1.red();
              texture[4*i+1] = (GLubyte) c1.green();
              texture[4*i+2] = (GLubyte) c1.blue();
              texture[4*i+3] = (GLubyte) 255;
            }
            // Bind the texture created above, set textures on. (From Red Book)

            glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

            glGenTextures(1, &texName);
            glBindTexture(GL_TEXTURE_1D, texName);

            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

            glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, static_cast<int>(n_colours), 0,
                         GL_RGBA, GL_UNSIGNED_BYTE, texture);

            glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

            glBindTexture(GL_TEXTURE_1D, texName);
            glEnable(GL_TEXTURE_1D);
            free(texture);
          }
        }

        if (s == ltsManager->selectedState())
        {
          c = GlUtil::Color::blend<GlUtil::Color::ColorMode::RGB>(c, SELECT_COLOR, SELECT_BLEND, 0);
        }

        glColor4ub(c.red(),c.green(),c.blue(),255);
        QVector3D p = s->getPositionAbs();

        glPushMatrix();
        glTranslatef(p.x(), p.y(), p.z());
        glScalef(ns, ns, ns);

        glPushName(s->getID());

        primitiveFactory.drawSimpleSphere();
        glPopName();
        glPopMatrix();

        if (markManager->isMarked(s))
        {
          glDisable(GL_TEXTURE_1D);

        }

        glDeleteTextures(1, &texName);
      }
    }
    Cluster* desc;
    for (int i = 0; i < root->getNumDescendants(); ++i)
    {
      desc = root->getDescendant(i);
      if (desc != NULL)
      {
        drawStates(desc, simulating);
      }
    }
  }
}

// ------------- TRANSITIONS ---------------------------------------------------

void Visualizer::drawTransitions(bool draw_fp,bool draw_bp)
{
  if (!ltsManager->lts())
  {
    return;
  }
  if (!draw_fp && !draw_bp)
  {
    return;
  }

  computeAbsPos();
  drawTransitions(ltsManager->lts()->getInitialState()->getCluster(),draw_fp,draw_bp);
}

void Visualizer::drawTransitions(Cluster* root,bool disp_fp,bool disp_bp)
{
  if (ltsManager->lts()->getZoomLevel() == root->getState(0)->getZoomLevel())
  {
    for (int j = 0; j < root->getNumStates(); ++j)
    {
      State* s = root->getState(j);
      // draw the outgoing transitions of s
      for (int i = 0; i < s->getNumOutTransitions(); ++i)
      {
        Transition* outTransition = s->getOutTransition(i);
        State* endState = outTransition->getEndState();

        if (ltsManager->lts()->getZoomLevel() == endState->getZoomLevel())
        {
          // Draw transition from root to endState
          if (disp_bp && outTransition->isBackpointer())
          {
            QColor from, to;
            if (markManager->isMarked(outTransition))
            {
              from = Settings::instance().markedColor.value();
              to = Settings::instance().markedColor.value();
            }
            else
            {
              from = Settings::instance().upEdgeColorFrom.value();
              to = Settings::instance().upEdgeColorTo.value();
            }
            // TODO: Use settings->longinterpolationUpEdge
            drawBackPointer(s,from,endState,to);
          } 
          else
          if (disp_fp && !outTransition->isBackpointer())
          {
            if (markManager->isMarked(outTransition))
            {
              QColor c = Settings::instance().markedColor.value();
              glColor4ub(c.red(),c.green(),c.blue(),255);
            }
            else
            {
              QColor c = Settings::instance().downEdgeColor.value();
              glColor4ub(c.red(),c.green(),c.blue(),255);
            }
            drawForwardPointer(s,endState);
          }
        }
      }
      // draw a loop if s has a loop
      if (disp_fp && s->getNumLoops() > 0)
      {
        bool hasMarkedLoop = false;
        for (int i = 0; i < s->getNumLoops() && !hasMarkedLoop; ++i)
        {
          if (markManager->isMarked(s->getLoop(i)))
          {
            hasMarkedLoop = true;
          }
        }
        if (hasMarkedLoop)
        {
          QColor c = Settings::instance().markedColor.value();
          glColor4ub(c.red(),c.green(),c.blue(),255);
        }
        else
        {
          QColor c = Settings::instance().downEdgeColor.value();
          glColor4ub(c.red(),c.green(),c.blue(),255);
        }
        drawLoop(s);
      }
    }
    for (int i = 0; i < root->getNumDescendants(); ++i)
    {
      Cluster* desc = root->getDescendant(i);
      if (desc != NULL)
      {
        drawTransitions(desc,disp_fp,disp_bp);
      }
    }
  }
}

void Visualizer::drawSimTransitions(bool draw_fp, bool draw_bp,
                                    QList<Transition*> transHis,
                                    QList<Transition*> posTrans,
                                    Transition* chosenTrans)
{
  computeAbsPos();

  // Draw the historical transitions.
  for (int i = 0; i < transHis.size(); ++i)
  {
    Transition* currTrans = transHis[i];
    State* beginState = currTrans->getBeginState();
    State* endState = currTrans->getEndState();

    if (ltsManager->lts()->getZoomLevel() == beginState->getZoomLevel() &&
        ltsManager->lts()->getZoomLevel() == endState->getZoomLevel())
    {
      // Draw transition from beginState to endState
      if (currTrans->isBackpointer() && draw_bp)
      {
        QColor color, black(0, 0, 0);
        if (markManager->isMarked(currTrans))
        {
          color = Settings::instance().markedColor.value();
        }
        else
        {
          color = Settings::instance().simPrevColor.value();
        }
        drawBackPointer(beginState, black, endState, color);
      }
      if (!currTrans->isBackpointer() && draw_fp)
      {
        if (markManager->isMarked(currTrans))
        {
          QColor c = Settings::instance().markedColor.value();
          glColor4ub(c.red(), c.green(), c.blue(), 255);
        }
        else
        {
          QColor transColor = Settings::instance().simPrevColor.value();
          glColor4ub(transColor.red(), transColor.green(), transColor.blue(), 255);
        }
        if (!currTrans->isSelfLoop())
        {
          drawForwardPointer(beginState, endState);
        }
        else
        {
          drawLoop(beginState);
        }
      }
    }
  }

  // Draw the possible transitions from the current state, as well as the state
  // they lead into
  for (int i = 0; i < posTrans.size(); ++i)
  {
    Transition* currTrans = posTrans[i];
    State* beginState = currTrans->getBeginState();
    State* endState = currTrans->getEndState();


    if (ltsManager->lts()->getZoomLevel() == beginState->getZoomLevel() &&
        ltsManager->lts()->getZoomLevel() == endState->getZoomLevel())
    {
      // Draw transition from beginState to endState
      if (currTrans->isBackpointer() && draw_bp)
      {
        QColor color, black(0, 0, 0);
        if (markManager->isMarked(currTrans))
        {
          color = Settings::instance().markedColor.value();
        }
        else
        if (currTrans == chosenTrans)
        {
          color = Settings::instance().simSelColor.value();
          glLineWidth(2.0);
        }
        else
        {
          color = Settings::instance().simPosColor.value();
        }
        drawBackPointer(beginState, black, endState, color);
        glLineWidth(1.0);
      } 
      else
      if (!currTrans->isBackpointer() && draw_fp)
      {
        if (currTrans == chosenTrans)
        {
          QColor c  = Settings::instance().simSelColor.value();
          glColor4ub(c.red(), c.green(), c.blue(), 255);
          glLineWidth(2.0);
        }
        else
        {
          QColor c = Settings::instance().simPosColor.value();
          glColor4ub(c.red(), c.green(), c.blue(), 255);
        }
        if (markManager->isMarked(currTrans))
        {
          QColor c = Settings::instance().markedColor.value();
          glColor4ub(c.red(), c.green(), c.blue(), 255);
        }

        if (currTrans->isSelfLoop())
        {
          drawLoop(beginState);
        }
        else
        {
          drawForwardPointer(beginState, endState);
        }

        glLineWidth(1.0);
      }
    }
  }
}

void Visualizer::drawForwardPointer(State* startState, State* endState)
{
  QVector3D startPoint = startState->getPositionAbs();
  QVector3D endPoint = endState->getPositionAbs();

  glBegin(GL_LINES);
  glVertex3f(startPoint.x(), startPoint.y(), startPoint.z());
  glVertex3f(endPoint.x(), endPoint.y(), endPoint.z());
  glEnd();
}

void Visualizer::drawBackPointer(State* startState, QColor& startColor, State* endState, QColor& endColor)
{
  QVector3D startPoint = startState->getPositionAbs();
  QVector3D startControl = startState->getOutgoingControl();
  QVector3D endControl = endState->getIncomingControl();
  QVector3D endPoint = endState->getPositionAbs();

  GLfloat ctrlPts [4][3] =
  {
    { startPoint.x(), startPoint.y(), startPoint.z() },
    { startControl.x(), startControl.y(), startControl.z() },
    { endControl.x(), endControl.y(), endControl.z() },
    { endPoint.x(), endPoint.y(), endPoint.z() }
  };

  if (startState->isCentered() && endState->isCentered())
  {
    ctrlPts[1][0] = startPoint.x() * 1.25;
    ctrlPts[2][0] = startControl.x();
  }

  float t,it,b0,b1,b2,b3,x,y,z;
  int N = Settings::instance().quality.value();
  QColor _col;
  glBegin(GL_LINE_STRIP);

  for (int k = 0; k < N; ++k)
  {
    t  = (float)k / (N-1);
    it = 1.0f - t;
    b3 =      t *  t *  t;
    b2 = 3 *  t *  t * it;
    b1 = 3 *  t * it * it;
    b0 =     it * it * it;

    x = b0 * ctrlPts[0][0] +
        b1 * ctrlPts[1][0] +
        b2 * ctrlPts[2][0] +
        b3 * ctrlPts[3][0];

    y = b0 * ctrlPts[0][1] +
        b1 * ctrlPts[1][1] +
        b2 * ctrlPts[2][1] +
        b3 * ctrlPts[3][1];

    z = b0 * ctrlPts[0][2] +
        b1 * ctrlPts[1][2] +
        b2 * ctrlPts[2][2] +
        b3 * ctrlPts[3][2];
    if (Settings::instance().longInterpolationUpEdge.value())
    {
        _col = GlUtil::Color::lerp<GlUtil::Color::InterpolateMode::LONG,
                                         GlUtil::Color::ColorMode::HSL>(
            startColor, endColor, t);
    }
    else
    {
      _col = GlUtil::Color::lerp<GlUtil::Color::InterpolateMode::SHORT,
                                 GlUtil::Color::ColorMode::HSL>(startColor,
                                                                endColor, t);
    }
    
    glColor4ub((GLubyte)_col.red(), (GLubyte)_col.green(), (GLubyte)_col.blue(), 255);
    glVertex3f(x, y, z);
  }

  glEnd();
}

void Visualizer::drawLoop(State* state)
{
  QVector3D statePoint = state->getPositionAbs();
  QVector3D startControl = state->getLoopControl1();
  QVector3D endControl = state->getLoopControl2();

  float t,it,b0,b1,b2,b3,x,y,z;
  int N = Settings::instance().quality.value();
  glBegin(GL_LINE_STRIP);
  for (int k = 0; k < N; ++k)
  {
    t  = (float)k / (N-1);
    it = 1.0f - t;
    b0 =      t *  t *  t;
    b1 = 3 *  t *  t * it;
    b2 = 3 *  t * it * it;
    b3 =     it * it * it;

    x = b0 * statePoint.x()  +
        b1 * startControl.x() +
        b2 * endControl.x() +
        b3 * statePoint.x();

    y = b0 * statePoint.y() +
        b1 * startControl.y() +
        b2 * endControl.y() +
        b3 * statePoint.y();

    z = b0 * statePoint.z() +
        b1 * startControl.z() +
        b2 * endControl.z() +
        b3 * statePoint.z();
    glVertex3f(x, y, z);
  }
  glEnd();
}

void Visualizer::exportToText(std::string filename)
{
  if (!ltsManager->lts())
  {
    return;
  }
  std::map< Cluster*, unsigned int > clus_id;
  unsigned int N = 0;
  float ch = Settings::instance().clusterHeight.value();
  Cluster_iterator ci(ltsManager->lts());
  std::ofstream file(filename.c_str());
  if (!file)
  {
    return ;
  }

  file << "Legend:\n";
  file << "  cone(x,y,z): a cone with top radius x, base radius y and height z\n";
  file << "  sphere(x):   a sphere with radius x\n";
  file << "  at angle(x): cluster is on the rim of the base of its parent cluster at angle x\n";
  file << "               If cluster is a cone, the center of its top side is placed there.\n";
  file << "               If cluster is a sphere, the center of the sphere is placed there.\n";
  file << "---------------------------------------------------------------------------------\n";

  for (; !ci.is_end(); ++ci)
  {
    file << "rank " << (*ci)->getRank() << ": cluster " << N << " ";
    if ((*ci)->hasDescendants())
    {
      file << "cone(" << (*ci)->getTopRadius() << "," << (*ci)->getBaseRadius()
           << "," << ch << ") ";
    }
    else
    {
      file << "sphere(" << (*ci)->getTopRadius() << ") ";
    }
    if ((*ci)->isCentered())
    {
      file << "centered";
    }
    else
    {
      file << "at angle(" << (*ci)->getPosition() << ")";
    }
    if ((*ci)->getAncestor() != NULL)
    {
      file << " below parent cluster " << clus_id[(*ci)->getAncestor()];
    }
    file << std::endl;
    clus_id[*ci] = N;
    ++N;
  }
}
