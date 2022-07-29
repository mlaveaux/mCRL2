// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "state.h"
#include "transition.h"


State::State(int aid):
  cluster(NULL),
  id(aid),
  zoomLevel(0),
  positionAngle(-1.0f),
  positionRadius(0.0f),
  rank(0),
  simulationCount(0)
{}

State::~State()
{
  unsigned int i;
  for (i = 0; i < outTransitions.size(); ++i)
  {
    delete outTransitions[i];
  }
  for (i = 0; i < loops.size(); ++i)
  {
    delete loops[i];
  }
}

void State::addInTransition(Transition* trans)
{
  inTransitions.push_back(trans);
}

void State::addOutTransition(Transition* trans)
{
  outTransitions.push_back(trans);
}

void State::addLoop(Transition* trans)
{
  loops.push_back(trans);
}

bool State::addMatchedRule(MarkRuleIndex index)
{
  matchedRules.insert(index);
  return matchedRules.size() == 1;
}

bool State::removeMatchedRule(MarkRuleIndex index)
{
  return matchedRules.erase(index) > 0;
}

bool State::isDeadlock() const
{
  return (outTransitions.size() + loops.size() == 0);
}

std::size_t State::getID()
{
  return id;
}

void State::setID(std::size_t i)
{
  id = i;
}

int State::getRank() const
{
  return rank;
}

void State::setRank(int r)
{
  rank = r;
}

bool State::isCentered() const
{
  return (positionAngle < -0.9f);
}

void State::center()
{
  positionRadius = 0.0f;
  positionAngle = -1.0f;
}

float State::getPositionAngle() const
{
  return positionAngle;
}

float State::getPositionRadius() const
{
  return positionRadius;
}

QVector3D State::getPositionAbs() const
{
  return positionAbs;
}

QVector3D State::getOutgoingControl() const
{
  return outgoingControl;
}

QVector3D State::getIncomingControl() const
{
  return incomingControl;
}

QVector3D State::getLoopControl1() const
{
  return loopControl1;
}

QVector3D State::getLoopControl2() const
{
  return loopControl2;
}

void State::setPositionRadius(float r)
{
  positionRadius = r;
}

void State::setPositionAngle(float a)
{
  positionAngle = a;
}

void State::setPositionAbs(const QVector3D& p)
{
  positionAbs = p;
}

void State::setOutgoingControl(const QVector3D& p)
{
  outgoingControl = p;
}

void State::setIncomingControl(const QVector3D& p)
{
  incomingControl = p;
}

void State::setLoopControl1(const QVector3D& p)
{
  loopControl1 = p;
}

void State::setLoopControl2(const QVector3D& p)
{
  loopControl2 = p;
}

Cluster* State::getCluster() const
{
  return cluster;
}

void State::setCluster(Cluster* c)
{
  cluster = c;
}

Transition* State::getInTransition(int i) const
{
  return inTransitions[i];
}

int State::getNumInTransitions() const
{
  return static_cast<int>(inTransitions.size());
}

Transition* State::getOutTransition(int i) const
{
  return outTransitions[i];
}

int State::getNumOutTransitions() const
{
  return static_cast<int>(outTransitions.size());
}

Transition* State::getLoop(int i) const
{
  return loops[i];
}

int State::getNumLoops() const
{
  return static_cast<int>(loops.size());
}

int State::getZoomLevel() const
{
  return zoomLevel;
}

void State::setZoomLevel(const int level)
{
  zoomLevel = level;
}
