// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef STATE_H
#define STATE_H

#include "markmanager.h"
#include "vectors.h"

class Transition;
class Cluster;

class State
{
  public:
    State(int aid);
    ~State();
    void addInTransition(Transition* trans);
    void addLoop(Transition* trans);
    void addOutTransition(Transition* trans);
    void center();

    Cluster* getCluster() const;
    std::size_t getID();
    Transition* getInTransition(int i) const;
    Transition* getOutTransition(int i) const;
    Transition* getLoop(int i) const;
    int getNumInTransitions() const;
    int getNumOutTransitions() const;
    int getNumLoops() const;

    float getPositionAngle() const;
    float getPositionRadius() const;
    QVector3D getPositionAbs() const;
    QVector3D getOutgoingControl() const;
    QVector3D getIncomingControl() const;
    QVector3D getLoopControl1() const;
    QVector3D getLoopControl2() const;
    int getRank() const;
    bool isSimulated() const { return simulationCount > 0; }
    bool isCentered() const;
    bool isDeadlock() const;
    bool hasTextures() const;

    bool addMatchedRule(MarkRuleIndex index);
    bool removeMatchedRule(MarkRuleIndex index);
    const std::set<MarkRuleIndex> &getMatchedRules() const { return matchedRules; }
    void clearMatchedRules() { matchedRules.clear(); }

    void setCluster(Cluster* c);
    void setPositionAngle(float a);
    void setPositionRadius(float r);
    void setPositionAbs(const QVector3D& p);
    void setOutgoingControl(const QVector3D& p);
    void setIncomingControl(const QVector3D& p);
    void setLoopControl1(const QVector3D& p);
    void setLoopControl2(const QVector3D& p);
    void setID(std::size_t i);
    void setRank(int r);
    void increaseSimulation() { simulationCount++; }
    void decreaseSimulation() { simulationCount--; }

    void setZoomLevel(const int level);
    int getZoomLevel() const;

  private:
    Cluster* cluster;
    std::size_t id;
    int zoomLevel;
    std::vector< Transition* > inTransitions;
    std::vector< Transition* > loops;
    std::vector< Transition* > outTransitions;
    std::set< MarkRuleIndex > matchedRules;
    float positionAngle;
    float positionRadius;
    QVector3D positionAbs;
    QVector3D outgoingControl;
    QVector3D incomingControl;
    QVector3D loopControl1;
    QVector3D loopControl2;
    int rank;
    int simulationCount;
};

#endif //STATE_H
