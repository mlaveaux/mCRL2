// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef INFODOCK_H
#define INFODOCK_H

#include "ui_infodock.h"
#include "markmanager.h"

class InfoDock : public QWidget
{
  Q_OBJECT

  public:
    InfoDock(QWidget *parent, LtsManager *ltsManager, MarkManager *markManager);

  protected slots:
    void ltsChanged();
    void markStatisticsChanged();
    void selectionChanged();

  private:
    Ui::InfoDock m_ui;
    LtsManager *m_ltsManager;
    MarkManager *m_markManager;
};

#endif
