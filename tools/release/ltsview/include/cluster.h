// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file cluster.h
/// \brief Header file of the Cluster class

#ifndef CLUSTER_H
#define CLUSTER_H

#include "markmanager.h"

class State;

// class for cluster comparison based on cluster radii
class Comp_BCRadius
{
  public:
    bool operator()(const Cluster*,const Cluster*) const;
};

// class for cluster comparison based on cluster volumes
class Comp_BCVolume
{
  public:
    bool operator()(const Cluster*,const Cluster*) const;
};

class Cluster
{
  public:
    Cluster(int r);
    ~Cluster();

    // Methods on descendants and ancestors.
    void addDescendant(Cluster* c);
    Cluster* getDescendant(int i) const;
    int getNumDescendants() const;
    bool hasDescendants() const;

    // Lets the cluster pretend that it does not have descendant i (for
    // zooming)
    void severDescendant(int i);

    // Restores the severance
    void healSeverance(int i);
    bool hasSeveredDescendants();
    void setAncestor(Cluster* c);
    Cluster* getAncestor() const;

    // Methods on states
    void addState(State* s);
    State* getState(int i) const;
    int getNumStates() const;
    int getNumDeadlocks() const;
    void addDeadlock();

    void addMatchedRule(MarkRuleIndex index);
    void removeMatchedRule(MarkRuleIndex index);
    const std::set<MarkRuleIndex> &getMatchedRules() const { return matchedRules; }
    int getNumMarkedStatesAll();
    int getNumMarkedStatesAny();
    void setNumMarkedStatesAll(int n);
    void setNumMarkedStatesAny(int n);
    int getNumMarkedTransitions();

    void positionStatesSpiral();

    // Methods on transitions
    void addActionLabel(int l);
    void resetActionMarks() { numMarkedTransitions = 0; }
    bool hasMarkedTransition() const;
    int setActionMark(int l,bool b);

    // General cluster information
    void center();
    void computeSizeAndPositions();
    void computeSizeAndPositions_FSM();
    float getBaseRadius() const;
    float getPosition() const;
    int getRank() const;
    int getPositionInRank() const;
    void setPositionInRank(int p);
    float getTopRadius() const;
    float getBCHeight() const;
    float getBCRadius() const;
    float getBCVolume() const;
    bool isCentered() const;
    void setPosition(float p);

    // Methods for visualization
    int getVisObject() const;
    void setVisObject(int vo);
    int getBranchVisObject(int i) const;
    int getNumBranchVisObjects() const;
    void addBranchVisObject(int vo);
    void clearBranchVisObjects();

    std::vector< Cluster* > descendants;
  private:
    std::map< int, int > actionLabelCounts;
    Cluster* ancestor;
    float baseRadius;
    int numDeadlocks;
    std::vector<std::vector<bool> > severedDescendants;
    unsigned int severedDescendantsC;
    int numMarkedStatesAll;
    int numMarkedStatesAny;
    int numMarkedTransitions;
    float position;
    int rank;
    int positionInRank;
    std::vector< State* > states;
    float topRadius;
    float bc_radius; // radius of the bounding cylinder that contains this cluster's subtree
    float bc_height; // height of the bounding cylinder that contains this cluster's subtree
    int visObject;
    std::vector< int > branchVisObjects;
    std::set< MarkRuleIndex > matchedRules;
};

#endif
