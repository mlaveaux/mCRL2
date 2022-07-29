// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "simdock.h"

SimDock::SimDock(QWidget *parent, LtsManager *ltsManager)
  : QWidget(parent),
    m_ltsManager(ltsManager)
{
  m_ui.setupUi(this);

  connect(m_ui.start, SIGNAL(clicked()), this, SLOT(start()));
  connect(m_ui.stop, SIGNAL(clicked()), this, SLOT(stop()));
  connect(m_ui.backtrace, SIGNAL(clicked()), this, SLOT(backtrace()));
  connect(m_ui.reset, SIGNAL(clicked()), this, SLOT(reset()));
  connect(m_ui.trigger, SIGNAL(clicked()), this, SLOT(trigger()));
  connect(m_ui.undo, SIGNAL(clicked()), this, SLOT(undo()));
  connect(m_ui.transitionTable, SIGNAL(itemSelectionChanged()), this, SLOT(select()));
  connect(m_ui.transitionTable, SIGNAL(cellDoubleClicked(int, int)), this, SLOT(trigger()));

  connect(m_ltsManager, SIGNAL(simulationChanged()), this, SLOT(changed()));

  changed();
}

static inline QTableWidgetItem *item()
{
  QTableWidgetItem *item = new QTableWidgetItem();
  item->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
  return item;
}

void SimDock::changed()
{
  const std::unique_ptr<Simulation>& simulation = m_ltsManager->simulation();
  if (!simulation)
  {
    m_ui.start->setEnabled(false);
    m_ui.stop->setEnabled(false);
    m_ui.backtrace->setEnabled(false);
    m_ui.reset->setEnabled(false);
    m_ui.undo->setEnabled(false);
    m_ui.transitionTable->setRowCount(0);
  }
  else if (!simulation->isStarted())
  {
    m_ui.start->setEnabled(true);
    m_ui.stop->setEnabled(false);
    m_ui.backtrace->setEnabled(false);
    m_ui.reset->setEnabled(false);
    m_ui.undo->setEnabled(false);
    m_ui.transitionTable->setRowCount(0);
  }
  else
  {
    m_ui.start->setEnabled(false);
    m_ui.stop->setEnabled(true);
    m_ui.backtrace->setEnabled(true);
    m_ui.reset->setEnabled(true);
    m_ui.undo->setEnabled(simulation->canUndo());

    QList<Transition *> transitions = simulation->availableTransitions();
    State *currentState = simulation->currentState();

    m_ui.transitionTable->setRowCount(0);
    m_ui.transitionTable->setRowCount(transitions.size());
    for (int i = 0; i < transitions.size(); i++)
    {
      m_ui.transitionTable->setItem(i, 0, item());
      m_ui.transitionTable->setItem(i, 1, item());
      m_ui.transitionTable->item(i, 0)->setText(QString::fromStdString(simulation->lts().getActionLabel(transitions[i]->getLabel())));

      State *destination = transitions[i]->getEndState();
      QStringList assignments;
      for (std::size_t j = 0; j < simulation->lts().getNumParameters(); j++)
      {
        if (simulation->lts().getStateParameterValue(currentState, j) != simulation->lts().getStateParameterValue(destination, j))
        {
          assignments += QString::fromStdString(simulation->lts().getParameterName(j)) + QString(":=") + QString::fromStdString(simulation->lts().getStateParameterValueStr(destination, j));
        }
      }
      m_ui.transitionTable->item(i, 1)->setText(assignments.join(", "));
    }
    selectionChanged();
  }
}

void SimDock::selectionChanged()
{
  if (m_ltsManager->simulation() && m_ltsManager->simulation()->isStarted())
  {
    m_ui.trigger->setEnabled(true);
    QList<Transition *> transitions = m_ltsManager->simulation()->availableTransitions();
    for (int i = 0; i < transitions.size(); i++)
    {
      if (transitions[i] == m_ltsManager->simulation()->currentTransition())
      {
        m_ui.transitionTable->setRangeSelected(QTableWidgetSelectionRange(i, 0, i, 1), true);
      }
    }
  }
  else
  {
    m_ui.trigger->setEnabled(false);
  }
}

void SimDock::start()
{
  m_ltsManager->simulation()->start();
}

void SimDock::stop()
{
  m_ltsManager->simulation()->stop();
}

void SimDock::backtrace()
{
  m_ltsManager->simulation()->traceback();
}

void SimDock::reset()
{
  m_ltsManager->simulation()->stop();
  m_ltsManager->simulation()->start();
}

void SimDock::trigger()
{
  Transition *transition = m_ltsManager->simulation()->currentTransition();
  if (!transition)
  {
    QList<Transition *> transitions = m_ltsManager->simulation()->availableTransitions();
    if (!transitions.isEmpty())
    {
      transition = transitions[0];
    }
  }
  if (transition)
  {
    m_ltsManager->simulation()->followTransition(transition);
  }
}

void SimDock::undo()
{
  m_ltsManager->simulation()->undo();
}

void SimDock::select()
{
  QList<QTableWidgetSelectionRange> selection = m_ui.transitionTable->selectedRanges();
  if (selection.size() > 0)
  {
    m_ltsManager->simulation()->selectTransition(m_ltsManager->simulation()->availableTransitions()[selection[0].topRow()]);
  }
}
