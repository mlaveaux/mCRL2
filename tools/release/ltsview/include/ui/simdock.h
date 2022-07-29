// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef SIMDOCK_H
#define SIMDOCK_H

#include "ui_simdock.h"

#include "state.h"
#include "transition.h"

class SimDock : public QWidget
{
  Q_OBJECT

  public:
    SimDock(QWidget *parent, LtsManager *ltsManager);

  public slots:
    void changed();
    void selectionChanged();

  protected slots:
    void start();
    void stop();
    void backtrace();
    void reset();
    void trigger();
    void undo();
    void select();

  private:
    Ui::SimDock m_ui;
    LtsManager *m_ltsManager;
};

#endif
