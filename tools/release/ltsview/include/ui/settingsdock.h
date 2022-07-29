// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef SETTINGSDOCK_H
#define SETTINGSDOCK_H

#include "ui_settingsdock.h"

#include "settings.h"
#include <QDialog>

class ComboboxHandler : public QObject
{
  Q_OBJECT

  public:
    ComboboxHandler(QComboBox *combobox, Settings::SettingBool &setting);

  protected slots:
    void stateChanged(int state);
    void setState();

  private:
    QComboBox *m_combobox;
    Settings::SettingBool *m_setting;
};

class SettingsDock : public QWidget
{
  Q_OBJECT

  public:
    SettingsDock(QWidget *parent);

  private slots:
    void stateSizeChanged(int value);
    void setStateSize(float value);
    void clusterHeightChanged(int value);
    void setClusterHeight(float value);

  private:
    void setupSpinbox(QSpinBox *spinbox, Settings::SettingInt &setting);

  private:
    Ui::SettingsDock m_ui;
};

#endif
