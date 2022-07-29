// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file cluster.cpp
/// \brief Implementation of the Cluster class

#include <cmath>
#include "cluster.h"
#include "mathutils.h"
#include "state.h"

using namespace MathUtils;

bool Comp_BCRadius::operator()(const Cluster* c1,const Cluster* c2) const
{
  // The clause c1 != c2 is necessary in the following expression. If not
  // present, this comparison operator may no longer be irreflexive, i.e. it is
  // no longer guaranteed that Comp_BCRadius(x,x) returns false. This is due to
  // optimizations by the gcc compiler to (arbitrarily) store a float either in
  // memory or in a register. Depending on the architecture, these may have
  // different precisions, resulting in rounding errors.
  return c1 != c2 && c1->getBCRadius() < c2->getBCRadius();
}

bool Comp_BCVolume::operator()(const Cluster* c1,const Cluster* c2) const
{
  // See the explanation above on the c1 != c2 clause.
  return c1 != c2 && c1->getBCVolume() < c2->getBCVolume();
}

Cluster::Cluster(int r)
{
  ancestor = NULL;
  position = 0.0f;
  baseRadius = 0.0f;
  topRadius = 0.0f;
  numDeadlocks = 0;
  numMarkedStatesAll = 0;
  numMarkedStatesAny = 0;
  visObject = -1;
  numMarkedTransitions = 0;
  rank = r;
  severedDescendantsC = 0;
  bc_radius = 0.0f;
  bc_height = 0.0f;
  positionInRank = 0;

  for (std::size_t i = 0; i < descendants.size(); ++i)
  {
    std::vector<bool> dummy;
    severedDescendants.push_back(dummy);
  }
}

Cluster::~Cluster()
{
  actionLabelCounts.clear();
  descendants.clear();
  severedDescendants.clear();
  states.clear();
}

void Cluster::addState(State* s)
{
  states.push_back(s);
}

void Cluster::addActionLabel(int l)
{
  if (actionLabelCounts.find(l) == actionLabelCounts.end())
  {
    actionLabelCounts[l] = 1;
  }
  else
  {
    ++actionLabelCounts[l];
  }
}

void Cluster::setAncestor(Cluster* c)
{
  ancestor = c;
}

void Cluster::setPosition(float p)
{
  position = p;
}

void Cluster::addDescendant(Cluster* c)
{
  descendants.push_back(c);
  std::vector<bool> dummy;
  severedDescendants.push_back(dummy);
}

State* Cluster::getState(int i) const
{
  return states[i];
}

int Cluster::getNumStates() const
{
  return static_cast<int>(states.size());
}

Cluster* Cluster::getAncestor() const
{
  return ancestor;
}

int Cluster::getNumDescendants() const
{
  return static_cast<int>(descendants.size());
}

Cluster* Cluster::getDescendant(int i) const
{
  if (severedDescendants[i].size() > 0)
  {
    return NULL;
  }
  else
  {
    return descendants[i];
  }
}

bool Cluster::hasDescendants() const
{
  return descendants.size() != 0;
}


void Cluster::severDescendant(int i)
{
  severedDescendants[i].push_back(true);
  ++severedDescendantsC;
}


void Cluster::healSeverance(int i)
{
  if (severedDescendants[i].size() > 0)
  {
    severedDescendants[i].pop_back();
    --severedDescendantsC;
  }
}

bool Cluster::hasSeveredDescendants()
{
  return severedDescendantsC > 0;
}


float Cluster::getBCRadius() const
{
  return bc_radius;
}

float Cluster::getBCHeight() const
{
  return bc_height;
}

float Cluster::getBCVolume() const
{
  return bc_radius * bc_radius * static_cast<float>(PI) * bc_height;
}

float Cluster::getPosition() const
{
  return position;
}

void Cluster::center()
{
  position = -1.0f;
}

bool Cluster::isCentered() const
{
  return position < -0.9f;
}

int Cluster::getRank() const
{
  return rank;
}

int Cluster::getPositionInRank() const
{
  return positionInRank;
}

void Cluster::setPositionInRank(int p)
{
  positionInRank = p;
}

void Cluster::positionStatesSpiral()
{
  states[0]->center();
  float d = 0.25f;
  float r = d;
  unsigned int i = 1;
  float a,da;
  while (i < states.size())
  {
    a = 0.0f;
    da = rad_to_deg(d/r);
    while (i < states.size() && a < 360.0f - da)
    {
      states[i]->setPositionRadius(r);
      states[i]->setPositionAngle(a);
      ++i;
      a += da;
    }
    r += d;
  }
}

void Cluster::computeSizeAndPositions_FSM()
{
  // This process is described in Frank van Ham's Master's thesis, p. 24
  // Recurse into the tree (depth first)
  for (unsigned int i = 0; i < descendants.size(); ++i)
  {
    descendants[i]->computeSizeAndPositions_FSM();
  }

  /* Compute the cluster radius r such that all states fit nicely into the
   * cluster. Suppose the cluster has N states and every state is visualized by
   * a circle with radius 0.1. Suppose we want the total area of all states
   * to fill a quarter of the circle that represents the cluster. This way,
   * there is (hopefully) enough room left for transitions and space between the
   * states.
   * So: N * pi * 0.1^2 = 0.25 * pi * r^2
   * Hence: r = sqrt( N * 0.04 )
   */
  topRadius = states.size()/(2*static_cast<float>(PI));

  if (descendants.size() == 0)
  {
    baseRadius = topRadius;
    bc_radius = topRadius;
    bc_height = 1.0f;
  }
  else if (descendants.size() == 1)
  {
    Cluster* desc = *descendants.begin();
    baseRadius = desc->getTopRadius();
    bc_radius = std::max(topRadius,desc->getBCRadius());
    bc_height = desc->getBCHeight() + 1.0f;
    desc->center();

  }
  else     // descendants.size() > 1
  {
    // sort descendants by size in ascending order
    sort(descendants.begin(),descendants.end(),Comp_BCRadius());

    // determine whether a unique smallest descendant exists
    Cluster* smallest = descendants[0];
    Cluster* nextSmallest = descendants[1];
    bool uniqueSmallest = ((nextSmallest->getBCRadius()-smallest->getBCRadius()) /
                           smallest->getBCRadius()) > 0.01f;

    // determine whether a unique largest descendant exists
    Cluster* largest = descendants[descendants.size()-1];
    Cluster* nextLargest = descendants[descendants.size()-2];
    bool uniqueLargest = ((nextLargest->getBCRadius()-largest->getBCRadius()) /
                          largest->getBCRadius()) < -0.01f;

    // invariant: descendants in range [x,y) have not been assigned a position
    int x = 0;
    int y = static_cast<int>(descendants.size());
    float bcr_center = 0.0f;  // BC radius of largest descendant in center
    float bcr_rim = largest->getBCRadius();  // BC radius of largest descendant on rim
    if (uniqueLargest)
    {
      // center the largest descendant
      largest->center();
      --y;
      bcr_center = largest->getBCRadius();
      bcr_rim = nextLargest->getBCRadius();
    }
    if (uniqueSmallest && (!uniqueLargest || !smallest->hasDescendants()))
    {
      // center the smallest descendant
      smallest->center();
      ++x;
      if (!uniqueLargest)
      {
        bcr_center = smallest->getBCRadius();
      }
      if (y-x == 0)
      {
        bcr_rim = 0.0f;
      }
    }
    if (y-x == 1)
    {
      ++y;
      bcr_rim = largest->getBCRadius();
      if (smallest->isCentered())
      {
        bcr_center = smallest->getBCRadius();
      }
    }

    // compute the radius of the base of the cylinder and the cluster's size
    float min_radius1 = 0.0f;
    if (y-x > 1)
    {
      min_radius1 = (float)(bcr_rim / sin(PI / (y-x)));
    }
    float min_radius2 = bcr_center;
    if (y-x > 0)
    {
      min_radius2 += bcr_rim + 0.01f;
    }
    baseRadius = std::max(min_radius1,min_radius2);
    bc_radius = std::max(min_radius1 + bcr_rim,min_radius2);
    bc_radius = std::max(topRadius,bc_radius);
    bc_height = 0.0f;
    for (unsigned int i = 0; i < descendants.size(); ++i)
    {
      if (descendants[i]->getBCHeight() > bc_height)
      {
        bc_height = descendants[i]->getBCHeight();
      }
    }
    bc_height += 1.0f;

    // Divide the remaining descendants over the rim of the circle.
    float angle = 360.0f / (y-x);
    int i = 0;
    while (x+i != y)
    {
      descendants[x+i]->setPosition(i*angle);
      ++i;
    }
  }
}

void Cluster::computeSizeAndPositions()
{
  // This process is described in Frank van Ham's Master's thesis, p. 24
  // Recurse into the tree (depth first)
  for (unsigned int i = 0; i < descendants.size(); ++i)
  {
    descendants[i]->computeSizeAndPositions();
  }

  /* Compute the cluster radius r such that all states fit nicely into the
   * cluster. Suppose the cluster has N states and every state is visualized by
   * a circle with radius 0.1. Suppose we want the total area of all states
   * to fill a quarter of the circle that represents the cluster. This way,
   * there is (hopefully) enough room left for transitions and space between the
   * states.
   * So: N * pi * 0.1^2 = 0.25 * pi * r^2
   * Hence: r = sqrt( N * 0.04 )
   */
  topRadius = sqrt(states.size()*0.04f);

  if (descendants.size() == 0)
  {
    baseRadius = topRadius;
    bc_radius = topRadius;
    bc_height = 1.0f;
  }
  else if (descendants.size() == 1)
  {
    Cluster* desc = *descendants.begin();
    baseRadius = desc->getTopRadius();
    bc_radius = std::max(topRadius,desc->getBCRadius());
    bc_height = desc->getBCHeight() + 1.0f;
    desc->center();
  }
  else     // descendants.size() > 1
  {
    // sort descendants by size in ascending order
    sort(descendants.begin(),descendants.end(),Comp_BCVolume());

    // determine whether a unique smallest descendant exists
    Cluster* smallest = descendants[0];
    Cluster* nextSmallest = descendants[1];

    bool uniqueSmallest = ((nextSmallest->getBCVolume()-smallest->getBCVolume()) /
                           smallest->getBCVolume()) > 0.01f;

    // determine whether a unique largest descendant exists
    Cluster* largest = descendants[descendants.size()-1];
    Cluster* nextLargest = descendants[descendants.size()-2];

    bool uniqueLargest = ((nextLargest->getBCVolume()-largest->getBCVolume()) /
                          largest->getBCVolume()) < -0.01f;

    // invariant: descendants in range [x,y) have not been assigned a position
    int x = 0;
    int y = static_cast<int>(descendants.size());

    float bcr_center = 0.0f;  // BC radius of largest descendant in center
    float bcr_rim = largest->getBCRadius();  // BC radius of largest descendant on rim
    if (uniqueLargest)
    {
      // center the largest descendant
      largest->center();
      --y;
      bcr_center = largest->getBCRadius();
      bcr_rim = nextLargest->getBCRadius();
    }
    else
    {
      if (uniqueSmallest)
      {
        // center the smallest descendant
        smallest->center();
        ++x;
        bcr_center = smallest->getBCRadius();
      }
    }

    // compute the radius of the base of the cylinder and the cluster's size
    float minRimRadius = 0.0f;
    if (y-x > 1)
    {
      minRimRadius = (float)(bcr_rim / sin(PI / (y-x)));
    }
    baseRadius = std::max(bcr_center + bcr_rim + 0.01f,minRimRadius);
    bc_radius = std::max(bcr_center + bcr_rim + 0.01f,minRimRadius +bcr_rim);
    bc_radius = std::max(topRadius,bc_radius);
    bc_height = 0.0f;
    for (unsigned int i = 0; i < descendants.size(); ++i)
    {
      if (descendants[i]->getBCHeight() > bc_height)
      {
        bc_height = descendants[i]->getBCHeight();
      }
    }
    bc_height += 1.0f;

    // Divide the remaining descendants over the rim of the circle. First take
    // the two largest unpositioned clusters and place them opposite to each
    // other, then do the same for the two smallest unpositioned clusters. Keep
    // repeating these steps until all clusters have been positioned. So if the
    // list of clusters (sorted ascending by size) is [ 0, 1, 2, 3, 4, 5 ] then
    // the clusters are placed in the following order (starting at angle 0 and
    // going counter-clockwise): 5, 0, 3, 4, 1, 2

    int i = 0;
    int h = (y-x) / 2 + (y-x) % 2;
    float angle = 360.0f / (y-x);
    while (x != y)
    {
      if (i % 2 == 1)
      {
        descendants[x]->setPosition(i*angle);
        ++x;
        if (x != y)
        {
          descendants[x]->setPosition((h+i)*angle);
          ++x;
        }
      }
      else
      {
        descendants[y-1]->setPosition(i*angle);
        --y;
        if (x != y)
        {
          descendants[y-1]->setPosition((h+i)*angle);
          --y;
        }
      }
      ++i;
    }
  }
}

float Cluster::getTopRadius() const
{
  return topRadius;
}

float Cluster::getBaseRadius() const
{
  return baseRadius;
}

bool Cluster::hasMarkedTransition() const
{
  return (numMarkedTransitions > 0);
}

void Cluster::addMatchedRule(MarkRuleIndex index)
{
  matchedRules.insert(index);
}

void Cluster::removeMatchedRule(MarkRuleIndex index)
{
  matchedRules.erase(index);
}

int Cluster::setActionMark(int l,bool b)
{
  std::map<int,int>::iterator li = actionLabelCounts.find(l);
  if (li != actionLabelCounts.end())
  {
    if (b)
    {
      numMarkedTransitions += li->second;
    }
    else
    {
      numMarkedTransitions -= li->second;
    }
  }
  return numMarkedTransitions;
}

int Cluster::getNumDeadlocks() const
{
  return numDeadlocks;
}

void Cluster::addDeadlock()
{
  ++numDeadlocks;
}

int Cluster::getVisObject() const
{
  return visObject;
}

void Cluster::setVisObject(int vo)
{
  visObject = vo;
}

int Cluster::getBranchVisObject(int i) const
{
  return branchVisObjects[i];
}

int Cluster::getNumBranchVisObjects() const
{
  return static_cast<int>(branchVisObjects.size());
}

void Cluster::addBranchVisObject(int vo)
{
  branchVisObjects.push_back(vo);
}

void Cluster::clearBranchVisObjects()
{
  branchVisObjects.clear();
}

int Cluster::getNumMarkedTransitions()
{
  return numMarkedTransitions;
}

int Cluster::getNumMarkedStatesAll()
{
  return numMarkedStatesAll;
}

int Cluster::getNumMarkedStatesAny()
{
  return numMarkedStatesAny;
}

void Cluster::setNumMarkedStatesAll(int n)
{
  numMarkedStatesAll = n;
}

void Cluster::setNumMarkedStatesAny(int n)
{
  numMarkedStatesAny = n;
}
