// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MARKSTATERULEDIALOG_H
#define MARKSTATERULEDIALOG_H

#include "ui_markstateruledialog.h"

#include <QColor>
#include <QDialog>

#include "lts.h"

class MarkStateRuleDialog : public QDialog
{
  Q_OBJECT

  public:
    MarkStateRuleDialog(QWidget *parent, LTS* lts, QColor color, int parameter, bool negated, const std::set<int>& values);
    QColor color() { return m_color; }
    int parameter() { return m_currentParameter; }
    bool negated() { return m_ui.relationList->item(1)->isSelected(); }
    std::set<int> values();

  protected slots:
    void colorClicked();
    void setColor(QColor color);
    void parameterSelected();

  private:
    Ui::MarkStateRuleDialog m_ui;
    LTS *m_lts;
    QColor m_color;
    int m_currentParameter;
};

#endif
