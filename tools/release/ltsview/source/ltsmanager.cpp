// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "mcrl2/utilities/logger.h"
#include "state.h"
#include "transition.h"
#include <QEventLoop>

void LtsManagerHelper::loadLts(QString filename)
{
  emit loadingLts();

  m_lts = 0;

  LTS *lts = new LTS();
  try
  {
    if (!lts->readFromFile(filename.toStdString()))
    {
      delete lts;
      mCRL2log(mcrl2::log::error) << "Error loading file: " << filename.toStdString();
      emit finished();
      return;
    }
  }
  catch (const mcrl2::runtime_error& e)
  {
    delete lts;
    mCRL2log(mcrl2::log::error) << e.what();
    emit finished();
    return;
  }

  clusterStates(lts);
}

void LtsManagerHelper::clusterStates(LTS *lts)
{
  bool cyclic = Settings::instance().stateRankStyleCyclic.value();
  emit rankingStates();
  lts->rankStates(cyclic);
  emit clusteringStates();
  lts->clusterStates(cyclic);
  emit computingClusterInfo();
  lts->computeClusterInfo();
  positionClusters(lts);
}

void LtsManagerHelper::positionClusters(LTS *lts)
{
  emit positioningClusters();
  lts->positionClusters(Settings::instance().fsmStyle.value());
  positionStates(lts);
}

void LtsManagerHelper::positionStates(LTS *lts)
{
  emit positioningStates();
  lts->positionStates(Settings::instance().statePosStyleMultiPass.value());
  m_lts = lts;
  emit finished();
}

LtsManager::LtsManager(QObject *parent, QThread *atermThread):
  QObject(parent),
  m_lts(0),
  m_simulation(nullptr),
  m_selectedState(0),
  m_selectedCluster(0)
{
  m_helper.moveToThread(atermThread);

  connect(&m_helper, SIGNAL(loadingLts()), this, SIGNAL(loadingLts()));
  connect(&m_helper, SIGNAL(rankingStates()), this, SIGNAL(rankingStates()));
  connect(&m_helper, SIGNAL(clusteringStates()), this, SIGNAL(clusteringStates()));
  connect(&m_helper, SIGNAL(computingClusterInfo()), this, SIGNAL(computingClusterInfo()));
  connect(&m_helper, SIGNAL(positioningClusters()), this, SIGNAL(positioningClusters()));
  connect(&m_helper, SIGNAL(positioningStates()), this, SIGNAL(positioningStates()));

  connect(&Settings::instance().stateRankStyleCyclic, SIGNAL(changed(bool)), this, SLOT(clusterStates()));
  connect(&Settings::instance().fsmStyle, SIGNAL(changed(bool)), this, SLOT(positionClusters()));
  connect(&Settings::instance().statePosStyleMultiPass, SIGNAL(changed(bool)), this, SLOT(positionStates()));
}

State *LtsManager::currentSimulationState() const
{
  if (!m_simulation)
  {
    return 0;
  }
  else
  {
    return m_simulation->currentState();
  }
}

Transition *LtsManager::currentSimulationTransition() const
{
  if (!m_simulation)
  {
    return 0;
  }
  else
  {
    return m_simulation->currentTransition();
  }
}

QList<Transition *> LtsManager::simulationAvailableTransitions() const
{
  if (!m_simulation)
  {
    return QList<Transition *>();
  }
  else
  {
    return m_simulation->availableTransitions();
  }
}

bool LtsManager::openLts(QString filename)
{
  emit startStructuring();

  QEventLoop loop;
  connect(&m_helper, SIGNAL(finished()), &loop, SLOT(quit()));
  QMetaObject::invokeMethod(&m_helper, "loadLts", Q_ARG(QString, filename));
  loop.exec();

  LTS *lts = m_helper.lts();
  if (!lts)
  {
    emit stopStructuring();
    emit errorLoadingLts();
    return false;
  }

  emit ltsStructured();

  delete m_lts;
  m_lts = lts;
  m_simulation = std::unique_ptr<Simulation>(new Simulation(this, *lts));

  connect(m_simulation.get(), SIGNAL(changed()), this, SLOT(updateSimulationHistory()));
  connect(m_simulation.get(), SIGNAL(changed()), this, SLOT(unselectNonsimilated()));
  connect(m_simulation.get(), SIGNAL(changed()), this, SIGNAL(simulationChanged()));
  connect(m_simulation.get(), SIGNAL(selectionChanged()), this, SIGNAL(selectionChanged()));

  emit ltsChanged(lts);
  emit clustersChanged();
  emit clusterPositionsChanged();
  emit statePositionsChanged();
  emit simulationChanged();

  emit stopStructuring();
  return true;
}

void LtsManager::zoomInBelow(Cluster *cluster)
{
  unselect();
  m_lts = m_lts->zoomIntoBelow(cluster);
  emit ltsZoomed(m_lts);
}

void LtsManager::zoomInAbove(Cluster *cluster)
{
  unselect();
  m_lts = m_lts->zoomIntoAbove(cluster);
  emit ltsZoomed(m_lts);
}

bool LtsManager::zoomOut()
{
  LTS *lts = m_lts;
  m_lts = lts->zoomOut();
  if (lts != m_lts)
  {
    delete lts;
    emit ltsZoomed(m_lts);
    return true;
  }
  else
  {
    return false;
  }
}

void LtsManager::unselect()
{
  if (m_lts && (m_selectedState || m_selectedCluster))
  {
    m_selectedState = 0;
    m_selectedCluster = 0;
    emit selectionChanged();
  }
}

void LtsManager::selectState(State *state)
{
  if (m_lts)
  {
    m_selectedState = state;
    if (m_simulation->isStarted())
    {
      Transition *transition = transitionTo(state);
      if (transition)
      {
        m_simulation->selectTransition(transition);
      }
    }
    else
    {
      m_simulation->setInitialState(state);
    }
    emit selectionChanged();
  }
}

void LtsManager::selectCluster(Cluster *cluster)
{
  if (m_lts)
  {
    m_selectedCluster = cluster;
    emit selectionChanged();
  }
}

void LtsManager::simulateState(State *state)
{
  if (m_simulation)
  {
    Transition *transition = transitionTo(state);
    if (transition)
    {
      m_simulation->followTransition(transition);
    }
  }
}

void LtsManager::clusterStates()
{
  if (m_lts)
  {
    while (zoomOut())
    {
      // do nothing
    }

    emit startStructuring();

    QEventLoop loop;
    connect(&m_helper, SIGNAL(finished()), &loop, SLOT(quit()));
    QMetaObject::invokeMethod(&m_helper, "clusterStates", Q_ARG(LTS *, m_lts));
    loop.exec();

    emit ltsStructured();

    emit clustersChanged();
    emit clusterPositionsChanged();
    emit statePositionsChanged();

    emit stopStructuring();
  }
}

void LtsManager::positionClusters()
{
  if (m_lts)
  {
    emit startStructuring();

    QEventLoop loop;
    connect(&m_helper, SIGNAL(finished()), &loop, SLOT(quit()));
    QMetaObject::invokeMethod(&m_helper, "positionClusters", Q_ARG(LTS *, m_lts));
    loop.exec();

    emit ltsStructured();

    emit clusterPositionsChanged();
    emit statePositionsChanged();

    emit stopStructuring();
  }
}

void LtsManager::positionStates()
{
  if (m_lts)
  {
    emit startStructuring();

    QEventLoop loop;
    connect(&m_helper, SIGNAL(finished()), &loop, SLOT(quit()));
    QMetaObject::invokeMethod(&m_helper, "positionStates", Q_ARG(LTS *, m_lts));
    loop.exec();

    emit ltsStructured();

    emit statePositionsChanged();

    emit stopStructuring();
  }
}

void LtsManager::updateSimulationHistory()
{
  m_simulationStateHistory.clear();
  m_simulationTransitionHistory.clear();
  if (m_simulation)
  {
    QList<Transition *> history = m_simulation->history();
    for (int i = 0; i < history.size(); i++)
    {
      // TODO: this should only select states that are part of the current zoomed LTS!
      //if (history[i]->getBeginState()->
      {
        m_simulationStateHistory += history[i]->getBeginState();
        m_simulationTransitionHistory += history[i];
      }
    }
  }
}

void LtsManager::unselectNonsimilated()
{
  if (m_simulation->isStarted() == (transitionTo(selectedState()) == 0))
  {
    //unselect();
  }
}

Transition *LtsManager::transitionTo(State *state)
{
  QList<Transition *> transitions = m_simulation->availableTransitions();
  for (int i = 0; i < transitions.size(); i++)
  {
    if (transitions[i]->getEndState() == state)
    {
      return transitions[i];
    }
  }
  return 0;
}
