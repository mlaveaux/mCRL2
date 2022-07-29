// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef SAVEPICTUREDIALOG_H
#define SAVEPICTUREDIALOG_H

#include "ui_savepicturedialog.h"
#include <QDialog>

#include "ltscanvas.h"

class SavePictureDialog : public QDialog
{
  Q_OBJECT

  public:
    SavePictureDialog(QWidget *parent, LtsCanvas *canvas, QString filename);

  signals:
    void statusMessage(QString message);

  protected slots:
    void widthChanged(int value);
    void heightChanged(int value);
    void save();

  private:
    Ui::SavePictureDialog m_ui;
    LtsCanvas *m_canvas;
    QString m_filename;
    float m_width;
    float m_height;
    bool m_inChange;
};

#endif
