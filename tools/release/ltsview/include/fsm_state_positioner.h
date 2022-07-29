// Author(s): Bas Ploeger
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef FSM_STATE_POSITIONER_H
#define FSM_STATE_POSITIONER_H
#include <set>
#include <vector>
#include "state_positioner.h"
#include "vectors.h"

class ClusterSlotInfo;
class LTS;
class State;

class FSMStatePositioner: public StatePositioner
{
  public:
    FSMStatePositioner(LTS* l);
    ~FSMStatePositioner();
    void positionStates();

  private:
    std::map< Cluster*, ClusterSlotInfo* > slot_info;
    std::vector< State* > unpositioned_states;

    bool allStatesCentered(std::vector< State* > &states);
    void assignStateToSlot(State* state, int ring, int slot);
    void bottomUpPass();
    void getPredecessors(State* state, std::vector< State* >& predecessors);
    void getSuccessors(State* state, std::vector< State* >& successors);
    void requestStatePosition(State* state, QVector2D& position);
    void resolveUnpositioned();
    QVector2D sumStateVectorsInMultipleClusters(std::vector< State* > &states,
        float rim_radius);
    QVector2D sumStateVectorsInSingleCluster(std::vector< State* > &states);
    void topDownPass();
};

class ClusterSlotInfo
{
  public:
    ClusterSlotInfo(Cluster* cluster);
    ~ClusterSlotInfo() {}
    void findFarthestFreeSlot(int& ring, int& slot);
    void findNearestSlot(QVector2D& position, int& ring, int& slot);
    void findNearestFreeSlot(int& ring, int& slot);
    int getNumRings();
    int getNumSlots(int ring);
    void getPolarCoordinates(int ring, int slot, float& angle, float& radius);
    void occupySlot(int ring, int slot);

  private:
    struct Slot
    {
      Slot(int r, int s): ring(r), slot(s) {}
      int ring, slot;
    };

    struct Slot_less
    {
      bool operator()(const Slot& c1, const Slot& c2) const
      {
        return (c1.ring == c2.ring) ? c1.slot < c2.slot : c1.ring < c2.ring;
      }
    };

    typedef std::set< Slot, Slot_less > SlotSet;

    static const float MIN_DELTA_RING;
    static const float MIN_DELTA_SLOT;

    float delta_ring;
    std::vector< int > num_slots;
    SlotSet occupied_slots;

    QVector2D getVector(int ring, int slot);
};

#endif
