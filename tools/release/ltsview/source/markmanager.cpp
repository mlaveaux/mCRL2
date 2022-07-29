// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file markmanager.cpp
/// \brief Mark manager source file

#include "state.h"
#include "cluster.h"
#include "transition.h"


MarkManager::MarkManager(QObject *parent, LtsManager *ltsManager):
  QObject(parent),
  m_ltsManager(ltsManager),
  m_lts(0),
  m_markStyle(NO_MARKS),
  m_clusterMatchStyle(MATCH_ANY),
  m_stateMatchStyle(MATCH_ANY),
  m_markedStatesAny(0),
  m_markedStatesAll(0),
  m_markedTransitions(0),
  m_activeMarkRules(0)
{
  connect(ltsManager, SIGNAL(ltsChanged(LTS *)), this, SLOT(setLts(LTS *)));
  connect(ltsManager, SIGNAL(ltsZoomed(LTS *)), this, SLOT(setRelatedLts(LTS *)));
  connect(ltsManager, SIGNAL(clustersChanged()), this, SLOT(flushClusters()));
}

QList<MarkRuleIndex> MarkManager::markRules()
{
  QList<MarkRuleIndex> output;
  for (MarkRuleIndex i = m_markRules.begin(); i != m_markRules.end(); i++)
  {
    output += i;
  }
  return output;
}

int MarkManager::markedStates() const
{
  if (markStyle() == MARK_DEADLOCKS)
  {
    return m_lts->getNumDeadlocks();
  }
  else if (markStyle() == MARK_STATES)
  {
    if (stateMatchStyle() == MATCH_ALL)
    {
      return m_markedStatesAll;
    }
    else
    {
      return m_markedStatesAny;
    }
  }
  else
  {
    return 0;
  }
}

bool MarkManager::isMarked(State *state)
{
  if (markStyle() == MARK_DEADLOCKS)
  {
    return state->isDeadlock();
  }
  else if (markStyle() == MARK_STATES)
  {
    if (stateMatchStyle() == MATCH_ALL)
    {
      return (state->getMatchedRules().size() == (std::size_t)m_activeMarkRules);
    }
    else
    {
      return (!state->getMatchedRules().empty());
    }
  }
  else
  {
    return false;
  }
}

bool MarkManager::isMarked(Cluster *cluster)
{
  if (markStyle() == MARK_DEADLOCKS)
  {
    return cluster->getNumDeadlocks() > 0;
  }
  else if (markStyle() == MARK_STATES)
  {
    int marked = (stateMatchStyle() == MATCH_ALL) ? cluster->getNumMarkedStatesAll() : cluster->getNumMarkedStatesAny();
    int required = (clusterMatchStyle() == MATCH_ALL) ? cluster->getNumStates() : 1;
    return marked >= required;
  }
  else if (markStyle() == MARK_TRANSITIONS)
  {
    return cluster->hasMarkedTransition();
  }
  else
  {
    return false;
  }
}

bool MarkManager::isMarked(Transition *transition)
{
  return markStyle() == MARK_TRANSITIONS && m_markedActions[transition->getLabel()];
}

QList<QColor> MarkManager::markColors(State *state)
{
  QList<QColor> output;
  for (std::set<MarkRuleIndex>::iterator i = state->getMatchedRules().begin(); i != state->getMatchedRules().end(); i++)
  {
    output += (*i)->color;
  }
  return output;
}

QList<QColor> MarkManager::markColors(Cluster *cluster)
{
  QList<QColor> output;
  for (std::set<MarkRuleIndex>::iterator i = cluster->getMatchedRules().begin(); i != cluster->getMatchedRules().end(); i++)
  {
    output += (*i)->color;
  }
  return output;
}

void MarkManager::setLts(LTS *lts)
{
  setMarkStyle(NO_MARKS);
  setClusterMatchStyle(MATCH_ANY);
  setStateMatchStyle(MATCH_ANY);

  // It seems that this function is not necessary, and causes crashes when opening a new file when a file is already opened.
  // cleanLts();
  m_lts = lts;

  m_markedActions.clear();
  m_markedActions.resize(m_lts->getNumActionLabels());

  m_markRules.clear();

  m_markedStatesAny = 0;
  m_markedStatesAll = 0;
  m_markedTransitions = 0;
  m_activeMarkRules = 0;

  for (Cluster_iterator i = m_lts->getClusterIterator(); !i.is_end(); i++)
  {
    int states = (*i)->getNumStates();
    (*i)->setNumMarkedStatesAll((*i)->getNumStates());
    m_markedStatesAll += states;
  }

  emit ltsChanged();
  emit statisticsChanged();
  emit marksChanged();
}

void MarkManager::setRelatedLts(LTS *lts)
{
  m_lts = lts;
  flushMarkedStateNumbers();
  emit statisticsChanged();
}

void MarkManager::setMarkStyle(MarkStyle style)
{
  if (m_markStyle != style)
  {
    m_markStyle = style;
    emit markStyleChanged(style);
    emit marksChanged();
  }
}

void MarkManager::setClusterMatchStyle(MatchStyle style)
{
  assert(style != MATCH_MULTI);
  if (m_clusterMatchStyle != style)
  {
    m_clusterMatchStyle = style;
    emit clusterMatchStyleChanged(style);
    if (markStyle() == MARK_STATES)
    {
      emit marksChanged();
    }
  }
}

void MarkManager::setStateMatchStyle(MatchStyle style)
{
  if (m_stateMatchStyle != style)
  {
    m_stateMatchStyle = style;
    emit stateMatchStyleChanged(style);
    if (markStyle() == MARK_STATES)
    {
      emit marksChanged();
    }
  }
}

MarkRuleIndex MarkManager::addMarkRule(MarkRule rule)
{
  bool changed = false;
  MarkRuleIndex index = m_markRules.insert(m_markRules.end(), rule);
  if (rule.active)
  {
    applyRule(index);
    changed = true;
  }
  emit markRuleAdded(index);
  if (changed)
  {
    emit statisticsChanged();
    if (markStyle() == MARK_STATES)
    {
      emit marksChanged();
    }
  }
  return index;
}

void MarkManager::setMarkRule(MarkRuleIndex index, MarkRule rule)
{
  if (*index != rule)
  {
    bool changed = false;
    if (index->active)
    {
      unapplyRule(index);
      changed = true;
    }
    *index = rule;
    if (index->active)
    {
      applyRule(index);
      changed = true;
    }
    emit markRuleChanged(index);
    if (changed)
    {
      emit statisticsChanged();
      if (markStyle() == MARK_STATES)
      {
        emit marksChanged();
      }
    }
  }
}

void MarkManager::removeMarkRule(MarkRuleIndex index)
{
  bool changed = false;
  if (index->active)
  {
    unapplyRule(index);
    changed = true;
  }
  emit markRuleRemoved(index);
  m_markRules.erase(index);
  if (changed)
  {
    emit statisticsChanged();
    if (markStyle() == MARK_STATES)
    {
      emit marksChanged();
    }
  }
}

void MarkManager::setActionMarked(int action, bool marked)
{
  if (m_markedActions[action] != marked)
  {
    m_markedActions[action] = marked;

    m_markedTransitions = 0;
    for (Cluster_iterator i = m_lts->getClusterIterator(); !i.is_end(); i++)
    {
      (*i)->setActionMark(action, marked);
      m_markedTransitions += (*i)->getNumMarkedTransitions();
    }

    emit actionMarked(action, marked);
    emit statisticsChanged();
    if (markStyle() == MARK_TRANSITIONS)
    {
      emit marksChanged();
    }
  }
}

void MarkManager::flushClusters()
{
  for (Cluster_iterator i = m_lts->getClusterIterator(); !i.is_end(); i++)
  {
    Cluster *cluster = *i;

    int anyAdded = 0;
    int allAdded = 0;
    for (int j = 0; j < cluster->getNumStates(); j++)
    {
      State *state = cluster->getState(j);
      if (!state->getMatchedRules().empty())
      {
        anyAdded++;
      }
      if (state->getMatchedRules().size() == (std::size_t)m_activeMarkRules)
      {
        allAdded++;
      }
      for (std::set<MarkRuleIndex>::iterator k = state->getMatchedRules().begin(); k != state->getMatchedRules().end(); k++)
      {
        cluster->addMatchedRule(*k);
      }
    }
    cluster->setNumMarkedStatesAny(anyAdded);
    cluster->setNumMarkedStatesAll(allAdded);

    cluster->resetActionMarks();
    for (std::size_t j = 0; j < m_markedActions.size(); j++)
    {
      if (m_markedActions[j])
      {
        cluster->setActionMark(j, true);
      }
    }
  }
}

void MarkManager::flushMarkedStateNumbers()
{
  m_markedStatesAll = 0;
  m_markedStatesAny = 0;
  for (Cluster_iterator i = m_lts->getClusterIterator(); !i.is_end(); i++)
  {
    m_markedStatesAll += (*i)->getNumMarkedStatesAll();
    m_markedStatesAny += (*i)->getNumMarkedStatesAny();
  }
}

void MarkManager::cleanLts()
{
  if (m_lts)
  {
    for (State_iterator i = m_lts->getStateIterator(); !i.is_end(); i++)
    {
      (*i)->clearMatchedRules();
    }
  }
}

void MarkManager::applyRule(MarkRuleIndex index)
{
  for (Cluster_iterator i = m_lts->getClusterIterator(); !i.is_end(); i++)
  {
    Cluster *cluster = *i;
    bool matched = false;
    int anyAdded = 0;
    int allRemoved = 0;
    for (int j = 0; j < cluster->getNumStates(); j++)
    {
      State *state = cluster->getState(j);
      if (matchesRule(state, *index))
      {
        if (state->addMatchedRule(index))
        {
          anyAdded++;
        }
        matched = true;
      }
      else
      {
        if (state->getMatchedRules().size() == (std::size_t)m_activeMarkRules)
        {
          allRemoved++;
        }
      }
    }
    if (matched)
    {
      cluster->addMatchedRule(index);
    }
    cluster->setNumMarkedStatesAny(cluster->getNumMarkedStatesAny() + anyAdded);
    cluster->setNumMarkedStatesAll(cluster->getNumMarkedStatesAll() - allRemoved);
    m_markedStatesAny += anyAdded;
    m_markedStatesAll -= allRemoved;
  }
  m_activeMarkRules++;
}

void MarkManager::unapplyRule(MarkRuleIndex index)
{
  m_activeMarkRules--;
  for (Cluster_iterator i = m_lts->getClusterIterator(); !i.is_end(); i++)
  {
    Cluster *cluster = *i;
    bool matched = false;
    int allAdded = 0;
    int anyRemoved = 0;
    for (int j = 0; j < cluster->getNumStates(); j++)
    {
      State *state = cluster->getState(j);
      if (state->removeMatchedRule(index))
      {
        if (state->getMatchedRules().empty())
        {
          anyRemoved++;
        }
        matched = true;
      }
      else
      {
        if (state->getMatchedRules().size() == (std::size_t)m_activeMarkRules)
        {
          allAdded++;
        }
      }
    }
    if (matched)
    {
      cluster->removeMatchedRule(index);
    }
    cluster->setNumMarkedStatesAny(cluster->getNumMarkedStatesAny() - anyRemoved);
    cluster->setNumMarkedStatesAll(cluster->getNumMarkedStatesAll() + allAdded);
    m_markedStatesAny -= anyRemoved;
    m_markedStatesAll += allAdded;
  }
}

bool MarkManager::matchesRule(State *state, const MarkRule& rule) const
{
  return rule.negated ^ (rule.values.count(m_lts->getStateParameterValue(state, rule.parameter))>0);
}
