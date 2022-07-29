// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef SIMULATION_H
#define SIMULATION_H

#include <QList>
#include <QObject>

class LTS;
class State;
class Transition;
class LtsManager;

class Simulation : public QObject
{
  Q_OBJECT

  public:
    Simulation(QObject *parent, LTS& lts);
    ~Simulation();
    void operator=(const Simulation &other);
    const LTS& lts() const { return m_ltsRef; }
    bool isStarted() const { return m_currentState != 0; }
    State *initialState() const { return m_initialState; }
    State *currentState() const { return m_currentState; }
    Transition *currentTransition() const { return m_currentTransition; }
    QList<Transition *> history() const { return m_history; }
    QList<Transition *> availableTransitions() const;
    bool canUndo() const { return !m_history.isEmpty(); }

  public slots:
    void start();
    void stop();
    void setInitialState(State* initialState) { m_initialState = initialState; }
    void selectTransition(Transition *transition);
    void followTransition(Transition *transition);
    void undo();
    void traceback();
    bool loadTrace(QString filename);

  signals:
    void started();
    void stopped();
    void changed();
    void selectionChanged();

  private:
    LTS& m_ltsRef;
    State *m_initialState;
    State *m_currentState;
    Transition *m_currentTransition;
    QList<Transition *> m_history;
};

#endif
