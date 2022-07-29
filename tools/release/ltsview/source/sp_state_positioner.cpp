// Author(s): Bas Ploeger
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include "cluster.h"
#include "sp_state_positioner.h"
#include "mathutils.h"
#include "rtree.h"
#include "state.h"
#include "transition.h"


class ClusterStatePositioner
{
  public:
    ClusterStatePositioner(Cluster* c);

    virtual ~ClusterStatePositioner()
    { }

    virtual void positionStates()
    { }

  protected:
    static const float MIN_DELTA_RING;
    static const float MIN_DELTA_SLOT;

    Cluster* cluster;
    float delta_ring;
    std::vector< int > num_ring_slots;
};

const float ClusterStatePositioner::MIN_DELTA_RING = 0.22f;
const float ClusterStatePositioner::MIN_DELTA_SLOT = 0.22f;

class NodeClusterStatePositioner: public ClusterStatePositioner
{
  public:
    NodeClusterStatePositioner(Cluster* c):
      ClusterStatePositioner(c), slot_rtree(NULL)
    { }

    ~NodeClusterStatePositioner();
    void positionStates();

  private:
    void assignStateToNearestSlot(State* state, const QVector2D& position);
    void buildRTree();
    QVector2D sumSuccessorStateVectors(State* state);

    RTree* slot_rtree;
};


class LeafClusterStatePositioner: public ClusterStatePositioner
{
  public:
    LeafClusterStatePositioner(Cluster* c):
      ClusterStatePositioner(c)
    { }

    ~LeafClusterStatePositioner()
    { }

    void positionStates();

  private:
    void computeNumRingStates();
    std::vector< int > num_ring_states;
};


ClusterStatePositioner::ClusterStatePositioner(Cluster* c):
  cluster(c)
{
  int num_rings = static_cast<int>(1 + std::floor(cluster->getTopRadius() / MIN_DELTA_RING));
  if (num_rings > 1)
  {
    delta_ring = cluster->getTopRadius() / static_cast<float>(num_rings - 1);
  }
  else
  {
    delta_ring = 1.0f;
  }
  num_ring_slots.assign(num_rings, 0);
  for (int ring = 0; ring < num_rings; ++ring)
  {
    int num_slots = static_cast<int>(std::floor(2.0f * MathUtils::PI * ring * delta_ring /
                                     MIN_DELTA_SLOT));
    num_ring_slots[ring] = (num_slots > 1) ? num_slots : 1;
  }
}

NodeClusterStatePositioner::~NodeClusterStatePositioner()
{
  if (slot_rtree != NULL)
  {
    delete slot_rtree;
  }
}

void NodeClusterStatePositioner::positionStates()
{
  buildRTree();
  for (int s = 0; s < cluster->getNumStates(); ++s)
  {
    State* state = cluster->getState(s);
    QVector2D position = sumSuccessorStateVectors(state);
    assignStateToNearestSlot(state, position);
  }
}

void NodeClusterStatePositioner::buildRTree()
{
  PackedRTreeBuilder rtree_builder = PackedRTreeBuilder();
  for (std::size_t ring = 0; ring < num_ring_slots.size(); ++ring)
  {
    float radius = ring * delta_ring;
    int num_slots = num_ring_slots[ring];
    float delta_deg = 360.0f / static_cast<float>(num_slots);
    for (int slot = 0; slot < num_slots; ++slot)
    {
      rtree_builder.addPoint(Vectors::fromPolar(delta_deg * slot, radius));
    }
  }
  rtree_builder.buildRTree();
  slot_rtree = rtree_builder.getRTree();
}

QVector2D NodeClusterStatePositioner::sumSuccessorStateVectors(State* state)
{
  QVector2D sum_vector = QVector2D(0, 0);
  for (int t = 0; t < state->getNumOutTransitions(); ++t)
  {
    State* successor = state->getOutTransition(t)->getEndState();
    if (successor->getRank() == state->getRank() + 1)
    {
      if (successor->getCluster()->isCentered())
      {
        if (!successor->isCentered())
        {
          sum_vector += Vectors::fromPolar(successor->getPositionAngle(),
                                            successor->getPositionRadius());
        }
      }
      else
      {
        sum_vector += Vectors::fromPolar(
                        successor->getCluster()->getPosition(), cluster->getBaseRadius());
        if (!successor->isCentered())
        {
          sum_vector += Vectors::fromPolar(successor->getPositionAngle() +
                                            successor->getCluster()->getPosition(),
                                            successor->getPositionRadius());
        }
      }
    }
  }
  return sum_vector;
}

void NodeClusterStatePositioner::assignStateToNearestSlot(State* state,
    const QVector2D& position)
{
  QVector2D assigned_position = position;
  slot_rtree->findNearestNeighbour(position);
  if (slot_rtree->hasFoundNeighbour())
  {
    assigned_position = slot_rtree->foundNeighbour();
  }
  float angle, radius;
  Vectors::toPolar(angle, radius, assigned_position);
  if (radius < 0.5f * delta_ring)
  {
    state->center();
  }
  else
  {
    state->setPositionAngle(angle);
    state->setPositionRadius(radius);
  }
  slot_rtree->deletePoint(assigned_position);
}

void LeafClusterStatePositioner::positionStates()
{
  computeNumRingStates();
  int state_begin = 0;
  if (num_ring_states[0] == 1)
  {
    cluster->getState(0)->center();
    ++state_begin;
  }
  for (std::size_t ring = 1; ring < num_ring_states.size(); ++ring)
  {
    float radius = delta_ring * ring;
    int num_states = num_ring_states[ring];
    float delta_deg = 360.0f / static_cast<float>(num_states);
    float angle = (ring % 2 == 1) ? 0.0f : 0.5f * delta_deg;
    for (int s = state_begin; s < state_begin + num_states; ++s)
    {
      cluster->getState(s)->setPositionAngle(angle);
      cluster->getState(s)->setPositionRadius(radius);
      angle += delta_deg;
    }
    state_begin += num_states;
  }
}

void LeafClusterStatePositioner::computeNumRingStates()
{
  int total_slots = 0;
  for (std::size_t ring = 0; ring < num_ring_slots.size(); ++ring)
  {
    total_slots += num_ring_slots[ring];
  }
  int todo_states = cluster->getNumStates();
  num_ring_states.assign(num_ring_slots.size(), 0);
  for (std::size_t ring = 0; ring < num_ring_slots.size(); ++ring)
  {
    if (ring == num_ring_slots.size() - 1)
    {
      num_ring_states[ring] = todo_states;
    }
    else
    {
      float rel_slots = static_cast<float>(num_ring_slots[ring]) /
                        static_cast<float>(total_slots);
      int num_states = MathUtils::round_to_int(cluster->getNumStates() *
                       rel_slots);
      num_ring_states[ring] = num_states;
      todo_states -= num_states;
    }
  }
}


SinglePassStatePositioner::SinglePassStatePositioner(LTS* l)
  : StatePositioner(l)
{
}

SinglePassStatePositioner::~SinglePassStatePositioner()
{
}

void SinglePassStatePositioner::positionStates()
{
  for (Reverse_cluster_iterator ci = lts->getReverseClusterIterator();
       !ci.is_end(); ++ci)
  {
    Cluster* cluster = *ci;
    ClusterStatePositioner* cs_positioner;
    if (cluster->getNumDescendants() == 0 || (cluster->getNumDescendants() == 1
        && cluster->getDescendant(0)->getNumStates() == 1))
    {
      cs_positioner = new LeafClusterStatePositioner(cluster);
    }
    else
    {
      cs_positioner = new NodeClusterStatePositioner(cluster);
    }
    cs_positioner->positionStates();
    delete cs_positioner;
  }
}

