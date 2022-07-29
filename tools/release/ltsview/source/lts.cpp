// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "mcrl2/lts/lts_algorithm.h"
#include "mcrl2/lts/lts_io.h"
#include "cluster.h"
#include "fsm_state_positioner.h"
#include "mathutils.h"
#include "sp_state_positioner.h"
#include "state.h"
#include "transition.h"

using namespace mcrl2::core;
using namespace mcrl2::lts;
using namespace MathUtils;

/**************************** Cluster iterators *******************************/

Cluster_iterator::Cluster_iterator(LTS* l)
{
  lts = l;
  rank = 0;
  cluster = 0;
  while (!is_end() && !is_valid())
  {
    next();
  }
}

void Cluster_iterator::operator++()
{
  if (!is_end())
  {
    next();
    while (!is_end() && !is_valid())
    {
      next();
    }
  }
}

Cluster* Cluster_iterator::operator*()
{
  return lts->clustersInRank[rank][cluster];
}

bool Cluster_iterator::is_valid()
{
  return (lts->clustersInRank[rank].size() > 0)
         && (lts->clustersInRank[rank][cluster] != NULL);
}

bool Cluster_iterator::is_end()
{
  return (rank == lts->getMaxRanks());
}

void Cluster_iterator::next()
{
  ++cluster;
  if (cluster >= lts->clustersInRank[rank].size())
  {
    cluster = 0;
    ++rank;
    while (rank < lts->getMaxRanks() &&
           lts->clustersInRank[rank].size() == 0)
    {
      ++rank;
    }
  }
}

Reverse_cluster_iterator::Reverse_cluster_iterator(LTS* l)
  : Cluster_iterator(l)
{
  rank = lts->getMaxRanks() - 1;
}

bool Reverse_cluster_iterator::is_end()
{
  return (rank == -1);
}

void Reverse_cluster_iterator::next()
{
  ++cluster;
  if (cluster >= lts->clustersInRank[rank].size())
  {
    cluster = 0;
    --rank;
    while (rank >= 0 && lts->clustersInRank[rank].size() == 0)
    {
      --rank;
    }
  }
}

/**************************** State_iterator **********************************/

State_iterator::State_iterator(LTS* l)
{
  lts = l;
  state_it = lts->states.begin();
}

State_iterator::~State_iterator()
{
}

void State_iterator::operator++()
{
  if (!is_end())
  {
    ++state_it;
  }
}

State* State_iterator::operator*()
{
  return *state_it;
}

bool State_iterator::is_end()
{
  return (state_it == lts->states.end());
}

/**************************** LTS *********************************************/

LTS::LTS()
{
  deadlockCount = -1;
  initialState = NULL;
  lastCluster = NULL;
  previousLevel = NULL;
  lastWasAbove = false;
  zoomLevel = 0;
}

LTS::LTS(LTS* parent, Cluster *target, bool fromAbove)
{
  // Copied/initial information
  lastWasAbove = fromAbove;
  previousLevel  = parent;
  zoomLevel = previousLevel->getZoomLevel() + 1;

  haveStateInfo = previousLevel->haveStateInfo;
  numActionLabels = previousLevel->numActionLabels;
  numParameters = previousLevel->numParameters;
  numStateLabels = previousLevel->numStateLabels;
  numStates = previousLevel->numStates;
  numTransitions = previousLevel->numTransitions;
  stateElementValues = previousLevel->stateElementValues;
  parameterNames = previousLevel->parameterNames;
  stateLabels = previousLevel->stateLabels;
  actionLabels = previousLevel->actionLabels;

  if (lastWasAbove)
  {
    initialState = previousLevel->getInitialState();

    lastCluster = target;

    Cluster* child = NULL;
    Cluster* parent = target;

    do
    {
      for (int i = 0; i < parent->getNumDescendants(); ++i)
      {
        if (child == NULL || parent->getDescendant(i) != child)
        {
          parent->severDescendant(i);
        }
      }

      addCluster(parent);
      child = parent;
      parent = child->getAncestor();
    }
    while (child != initialState->getCluster());

  }
  else
  {
    initialState = target->getState(0);
    addClusterAndBelow(target);
  }
}

LTS::~LTS()
{
  if (previousLevel == NULL)
  {
    // This LTS is the top level LTS, so delete all its contents.
    std::size_t i,r;
    std::vector< State* >::iterator li;
    for (li = states.begin(); li != states.end(); ++li)
    {
      delete *li;
    }
    states.clear();
    initialState = NULL;

    for (r = 0; r < clustersInRank.size(); ++r)
    {
      for (i = 0; i < clustersInRank[r].size(); ++i)
      {
        delete clustersInRank[r][i];
      }
    }
    clustersInRank.clear();
  }
  else
  {
    states.clear();
    clustersInRank.clear();
  }
}

Cluster_iterator LTS::getClusterIterator()
{
  return Cluster_iterator(this);
}

Reverse_cluster_iterator LTS::getReverseClusterIterator()
{
  return Reverse_cluster_iterator(this);
}

State_iterator LTS::getStateIterator()
{
  return State_iterator(this);
}

const std::string& LTS::getParameterName(std::size_t parindex) const
{
  assert(parindex < parameterNames.size());
  return parameterNames[parindex];
}

std::size_t LTS::getStateParameterValue(State* state,std::size_t param) const
{
  if (state->getID()<stateLabels.size() && param < stateLabels[state->getID()].size())
  { 
    return stateLabels[state->getID()][param];
  }
  return std::size_t(-1);
}

const std::string LTS::getStateParameterValueStr(State* state, std::size_t param) const
{
  std::size_t value = getStateParameterValue(state, param);
  if (value < stateElementValues[param].size())
  {
    return stateElementValues[param][value];
  }
  return std::string();
}

std::set<std::string> LTS::getClusterParameterValues(Cluster* cluster, std::size_t param) const
{
  std::set<std::string> result;
  for (int i = 0; i < cluster->getNumStates(); ++i)
  {
    result.insert(getStateParameterValueStr(cluster->getState(i),param));
  }
  return result;
}

bool LTS::readFromFile(const std::string& filename)
{
  assert(!previousLevel);

  mcrl2::lts::lts_fsm_t mcrl2_lts{};

  load_lts_as_fsm_file(filename, mcrl2_lts);

  // remove unreachable states
  reachability_check(mcrl2_lts, true);

  haveStateInfo = mcrl2_lts.has_state_info();
  numActionLabels = mcrl2_lts.num_action_labels();
  numParameters = mcrl2_lts.process_parameters().size();
  numStateLabels = mcrl2_lts.num_state_labels();
  numStates = mcrl2_lts.num_states();
  numTransitions = mcrl2_lts.num_transitions();

  stateElementValues.clear();
  stateElementValues.reserve(numParameters);
  for (std::size_t i = 0; i < numParameters; ++i)
  {
    stateElementValues.emplace_back(mcrl2_lts.state_element_values(i));
  }

  parameterNames.clear();
  parameterNames.reserve(numParameters);
  for (std::size_t i = 0; i < numParameters; ++i)
  {
    // in an .fsm file a parameter is a pair of strings.
    parameterNames.emplace_back(mcrl2_lts.process_parameter(i).first);
  }

  stateLabels.clear();
  stateLabels.reserve(numStateLabels);
  for (std::size_t i = 0; i < numStateLabels; ++i)
  {
    stateLabels.emplace_back(mcrl2_lts.state_label(i));
  }

  actionLabels.clear();
  actionLabels.reserve(numActionLabels);
  for (std::size_t i = 0; i < numActionLabels; ++i)
  {
    actionLabels.emplace_back(mcrl2_lts.action_label(i));
  }

  states.clear();
  states.reserve(mcrl2_lts.num_states());
  for (std::size_t i = 0; i < mcrl2_lts.num_states(); ++i)
  {
    states.emplace_back(new State(static_cast<int>(i)));
  }

  initialState = states[mcrl2_lts.initial_state()];

  const std::vector<transition> &trans = mcrl2_lts.get_transitions();
  for (const transition& r : trans)
  {
    State* s1 = states[r.from()];
    State* s2 = states[r.to()];
    Transition* t = new Transition(s1,s2,static_cast<int>(r.label()));
    if (s1 != s2)
    {
      s1->addOutTransition(t);
      s2->addInTransition(t);
    }
    else
    {
      s1->addLoop(t);
    }
  }
  return true;
}

int LTS::getNumActionLabels() const
{
  return static_cast<int>(numActionLabels);
}

std::size_t LTS::getNumParameters() const
{
  return numParameters;
}

const std::string& LTS::getActionLabel(int labindex) const
{
  return actionLabels[labindex];
}

void LTS::addCluster(Cluster* cluster)
{
  std::size_t rank = cluster->getRank();
  std::size_t pos = cluster->getPositionInRank();

  // Check to see if there is already a rank for this cluster
  if (clustersInRank.size() <= rank)
  {
    clustersInRank.resize(rank + 1);
  }

  if (clustersInRank[rank].size() <= pos)
  {
    clustersInRank[rank].resize(pos + 1);
  }

  clustersInRank[rank][pos] = cluster;

  for (int i = 0; i < cluster->getNumStates(); ++i)
  {
    State* state = cluster->getState(i);
    std::size_t state_id = state->getID();
    if (states.size() <= state_id)
    {
      states.resize(state_id + 1);
    }

    states[state_id] = state;
    state->setZoomLevel(zoomLevel);
  }
}

void LTS::addClusterAndBelow(Cluster* cluster)
{
  if (cluster != NULL)
  {
    addCluster(cluster);

    for (int i = 0; i < cluster->getNumDescendants(); ++i)
    {
      addClusterAndBelow(cluster->getDescendant(i));
    }
  }
}

State* LTS::getInitialState() const
{
  return initialState;
}

int LTS::getNumRanks() const
{
  int offset = 0;
  while (clustersInRank[offset].empty())
  {
    ++offset;
  }
  return static_cast<int>(clustersInRank.size()) - offset;
}

int LTS::getMaxRanks() const
{
  return static_cast<int>(clustersInRank.size());
}

int LTS::getNumClusters() const
{
  int result = 0;
  for (std::size_t i = 0; i < clustersInRank.size(); ++i)
  {
    result += static_cast<int>(clustersInRank[i].size());
  }
  return result;
}

int LTS::getNumDeadlocks()
{
  if (deadlockCount == -1)
  {
    // a value of -1 indicates that we have to compute it
    deadlockCount = 0;
    std::vector< State* >::iterator state_it;
    for (state_it  = states.begin(); state_it != states.end(); ++state_it)
    {
      if ((**state_it).isDeadlock())
      {
        ++deadlockCount;
      }
    }
  }
  return deadlockCount;
}

int LTS::getNumStates() const
{
  return static_cast<int>(numStates);
}

int LTS::getNumTransitions() const
{
  return static_cast<int>(numTransitions);
}

void LTS::clearRanksAndClusters()
{
  std::vector< State* >::iterator it;
  for (it = states.begin(); it != states.end(); ++it)
  {
    (*it)->setRank(-1);
    (*it)->setCluster(NULL);
  }

  for (Cluster_iterator ci = getClusterIterator(); !ci.is_end(); ++ci)
  {
    delete *ci;
  }

  std::vector< std::vector< Cluster* > > temp2;
  clustersInRank.swap(temp2);
}

void LTS::rankStates(bool cyclic)
{
  clearRanksAndClusters();
  int rankNumber = 0;

  std::vector< State* > nextRank,currRank;
  currRank.push_back(initialState);
  initialState->setRank(rankNumber);

  int i;
  std::vector< State* >::iterator it;
  State* s,*t;
  while (currRank.size() > 0)
  {
    nextRank.clear();
    // iterate over the states in this rank
    for (it  = currRank.begin(); it != currRank.end(); ++it)
    {
      s = *it;
      if (cyclic)
      {
        // iterate over all in-transitions of cs
        for (i = 0; i < s->getNumInTransitions(); ++i)
        {
          t = s->getInTransition(i)->getBeginState();
          if (t->getRank() == -1)
          {
            t->setRank(rankNumber+1);
            nextRank.push_back(t);
          }
        }
      }
      // iterate over all out-transitions of cs
      for (i = 0; i < s->getNumOutTransitions(); ++i)
      {
        t = s->getOutTransition(i)->getEndState();
        if (t->getRank() == -1)
        {
          t->setRank(rankNumber+1);
          nextRank.push_back(t);
        }
      }
    }
    currRank.swap(nextRank);
    ++rankNumber;
  }
}

void LTS::clusterStates(bool cyclic)
{
  Cluster* d = new Cluster(0);
  std::vector< Cluster* > cs;
  cs.push_back(d);
  d->setPositionInRank(0);
  clustersInRank.push_back(cs);
  clusterTree(initialState,d,cyclic);
}

void LTS::clusterTree(State* v,Cluster* c,bool cyclic)
{
  int h,i,j,r;
  State* w, *y;
  c->addState(v);
  v->setCluster(c);
  for (i = 0; i < v->getNumOutTransitions(); ++i)
  {
    w = v->getOutTransition(i)->getEndState();
    if (w->getCluster() == NULL && w->getRank() == v->getRank())
    {
      clusterTree(w,c,cyclic);
    }
  }
  for (i = 0; i < v->getNumInTransitions(); ++i)
  {
    w = v->getInTransition(i)->getBeginState();
    if (w->getCluster() == NULL && w->getRank() == v->getRank())
    {
      clusterTree(w,c,cyclic);
    }
  }
  if (cyclic)
  {
    for (i = 0; i < v->getNumInTransitions(); ++i)
    {
      w = v->getInTransition(i)->getBeginState();
      r = w->getRank();
      if (w->getCluster() == NULL && r == v->getRank()+1)
      {
        Cluster* d = new Cluster(r);
        if ((std::size_t)(r) >= clustersInRank.size())
        {
          std::vector< Cluster* > cs;
          clustersInRank.push_back(cs);
        }
        d->setPositionInRank(static_cast<int>(clustersInRank[r].size()));
        clustersInRank[r].push_back(d);
        d->setAncestor(c);
        c->addDescendant(d);
        clusterTree(w,d,cyclic);

        for (h = 0; h < d->getNumStates(); ++h)
        {
          y = d->getState(h);
          for (j = 0; j < y->getNumOutTransitions(); ++j)
          {
            w = y->getOutTransition(j)->getEndState();
            if (w->getCluster() == NULL && w->getRank() == v->getRank())
            {
              clusterTree(w,c,cyclic);
            }
          }
          for (j = 0; j < y->getNumInTransitions(); ++j)
          {
            w = y->getInTransition(j)->getBeginState();
            if (w->getCluster() == NULL && w->getRank() == v->getRank())
            {
              clusterTree(w,c,cyclic);
            }
          }
        }
      }
    }
  }
  for (i = 0; i < v->getNumOutTransitions(); ++i)
  {
    w = v->getOutTransition(i)->getEndState();
    r = w->getRank();
    if (w->getCluster() == NULL && r == v->getRank()+1)
    {
      Cluster* d = new Cluster(r);
      if ((std::size_t)(r) >= clustersInRank.size())
      {
        std::vector< Cluster* > cs;
        clustersInRank.push_back(cs);
      }
      d->setPositionInRank(static_cast<int>(clustersInRank[r].size()));
      clustersInRank[r].push_back(d);
      d->setAncestor(c);
      c->addDescendant(d);
      clusterTree(w,d,cyclic);

      for (h = 0; h < d->getNumStates(); ++h)
      {
        y = d->getState(h);
        if (cyclic)
        {
          for (j = 0; j < y->getNumOutTransitions(); ++j)
          {
            w = y->getOutTransition(j)->getEndState();
            if (w->getCluster() == NULL && w->getRank() == v->getRank())
            {
              clusterTree(w,c,cyclic);
            }
          }
        }
        for (j = 0; j < y->getNumInTransitions(); ++j)
        {
          w = y->getInTransition(j)->getBeginState();
          if (w->getCluster() == NULL && w->getRank() == v->getRank())
          {
            clusterTree(w,c,cyclic);
          }
        }
      }
    }
  }
}

void LTS::computeClusterInfo()
{
  State* s;
  Cluster* c;
  std::vector< State* >::iterator li;
  int t;
  for (li = states.begin(); li != states.end(); ++li)
  {
    s = *li;
    c = s->getCluster();
    if (s->isDeadlock())
    {
      c->addDeadlock();
    }
    for (t = 0; t < s->getNumOutTransitions(); ++t)
    {
      c->addActionLabel(s->getOutTransition(t)->getLabel());
    }
    for (t = 0; t < s->getNumLoops(); ++t)
    {
      c->addActionLabel(s->getLoop(t)->getLabel());
    }
  }
}

void LTS::positionClusters(bool fsmstyle)
{
  if (fsmstyle)
  {
    initialState->getCluster()->computeSizeAndPositions_FSM();
  }
  else
  {
    initialState->getCluster()->computeSizeAndPositions();
  }
  // position the initial state's cluster
  initialState->getCluster()->center();
}

void LTS::clearStatePositions()
{
  for (State_iterator state_it = getStateIterator(); !state_it.is_end();
       ++state_it)
  {
    (**state_it).center();
  }
}

void LTS::positionStates(bool multiPass)
{
  StatePositioner* state_positioner = NULL;
  if (multiPass)
  {
    state_positioner = new FSMStatePositioner(this);
  }
  else
  {
    state_positioner = new SinglePassStatePositioner(this);
  }
  state_positioner->positionStates();
  delete state_positioner;
}

LTS* LTS::zoomIntoAbove(Cluster *target)
{
  return new LTS(this, target, true);
}

LTS* LTS::zoomIntoBelow(Cluster *target)
{
  return new LTS(this, target, false);
}

LTS* LTS::zoomOut()
{
  if (previousLevel != NULL)
  {
    if (lastWasAbove)
    {
      Cluster* child = NULL;
      Cluster* parent = lastCluster;
      do
      {
        for (int i = 0; i < parent->getNumDescendants(); ++i)
        {
          parent->healSeverance(i);
        }
        child = parent;
        parent = child->getAncestor();
      }
      while (child != initialState->getCluster());
    }
    std::vector< State* >::iterator li;
    for (li = states.begin(); li != states.end(); ++li)
    {
      if (*li)
      {
        (*li)->setZoomLevel(zoomLevel - 1);
      }
    }
    return previousLevel;
  }
  else
  {
    return this;
  }
}

void LTS::setLastCluster(Cluster* c)
{
  lastCluster = c;
}

void LTS::setPreviousLevel(LTS* prev)
{
  previousLevel = prev;
}

LTS* LTS::getPreviousLevel() const
{
  return previousLevel;
}

void LTS::fromAbove()
{
  lastWasAbove = true;
}

int LTS::getZoomLevel() const
{
  return zoomLevel;
}

void LTS::setZoomLevel(const int level)
{
  zoomLevel = level;
}
