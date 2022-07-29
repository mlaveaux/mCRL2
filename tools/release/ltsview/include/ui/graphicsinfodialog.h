// Author(s): Ruben Vink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef GRAPHICSINFODIALOG_H
#define GRAPHICSINFODIALOG_H

#include "ui_graphicsinfodialog.h"

#include <QDialog>
#include <QLabel>
#include <QTextEdit>


class GraphicsInfoDialog : public QDialog
{
  Q_OBJECT

  public:
    GraphicsInfoDialog(QWidget *parent);

    QLabel *lbl_info;
    QTextEdit *txt_info;
  private:
    Ui::graphicsDebugInfoDialog m_ui;
};

#endif
