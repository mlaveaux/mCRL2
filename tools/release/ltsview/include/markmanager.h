// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file markmanager.h
/// \brief Mark manager header file

#ifndef MARKMANAGER_H
#define MARKMANAGER_H

#include "ltsmanager.h"

class Cluster;
class LTS;
struct MarkRule;
class State;
class Transition;

enum MarkStyle
{
  NO_MARKS,
  MARK_DEADLOCKS,
  MARK_STATES,
  MARK_TRANSITIONS
};

enum MatchStyle
{
  MATCH_ANY,
  MATCH_ALL,
  MATCH_MULTI
};

struct MarkRule
{
  bool active;
  QColor color;
  int parameter;
  bool negated;
  std::set<int> values;

  bool operator==(const MarkRule &other)
  {
    return active == other.active &&
           color == other.color &&
           parameter == other.parameter &&
           negated == other.negated &&
           values == other.values;
  }
  bool operator!=(const MarkRule &other) { return !(*this == other); }
};

typedef std::list<MarkRule>::iterator MarkRuleIndex;
inline bool operator<(const MarkRuleIndex &index1, const MarkRuleIndex &index2) { return &*index1 < &*index2; }

class MarkManager : public QObject
{
  Q_OBJECT

  public:
    MarkManager(QObject *parent, LtsManager *ltsManager);

    LTS *lts() const { return m_lts; }
    MarkStyle markStyle() const { return m_markStyle; }
    MatchStyle clusterMatchStyle() const { return m_clusterMatchStyle; }
    MatchStyle stateMatchStyle() const { return m_stateMatchStyle; }
    QList<MarkRuleIndex> markRules();
    MarkRule markRule(MarkRuleIndex index) const { return *index; }
    std::vector<bool> markedActions() const { return m_markedActions; }
    bool isActionMarked(int action) const { return m_markedActions[action]; }
    int markedStates() const;
    int markedTransitions() const { return markStyle() == MARK_TRANSITIONS ? m_markedTransitions : 0; }

    bool isMarked(State *state);
    bool isMarked(Cluster *cluster);
    bool isMarked(Transition *transition);

    QList<QColor> markColors(State *state);
    QList<QColor> markColors(Cluster *cluster);

  public slots:
    void setLts(LTS *lts);
    void setRelatedLts(LTS *lts);
    void setMarkStyle(MarkStyle style);
    void setClusterMatchStyle(MatchStyle style);
    void setStateMatchStyle(MatchStyle style);
    MarkRuleIndex addMarkRule(MarkRule rule);
    void setMarkRule(MarkRuleIndex index, MarkRule rule);
    void removeMarkRule(MarkRuleIndex index);
    void setActionMarked(int action, bool marked);
    void flushClusters();

  protected slots:
    void flushMarkedStateNumbers();
    void cleanLts();
    void applyRule(MarkRuleIndex index);
    void unapplyRule(MarkRuleIndex index);

  protected:
    bool matchesRule(State *state, const MarkRule &rule) const;

  signals:
    void ltsChanged();
    void markStyleChanged(MarkStyle style);
    void clusterMatchStyleChanged(MatchStyle style);
    void stateMatchStyleChanged(MatchStyle style);
    void markRuleAdded(MarkRuleIndex index);
    void markRuleChanged(MarkRuleIndex index);
    void markRuleRemoved(MarkRuleIndex index);
    void actionMarked(int action, bool marked);
    void statisticsChanged();
    void marksChanged();

  private:
    LtsManager *m_ltsManager;
    LTS* m_lts;
    MarkStyle m_markStyle;
    MatchStyle m_clusterMatchStyle;
    MatchStyle m_stateMatchStyle;
    std::list<MarkRule> m_markRules;
    std::vector<bool> m_markedActions;

    int m_markedStatesAny;
    int m_markedStatesAll;
    int m_markedTransitions;
    int m_activeMarkRules;
};

#endif
