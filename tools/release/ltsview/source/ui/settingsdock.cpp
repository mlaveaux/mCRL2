// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "ui/settingsdock.h"

ComboboxHandler::ComboboxHandler(QComboBox *combobox, Settings::SettingBool &setting):
  QObject(combobox),
  m_combobox(combobox),
  m_setting(&setting)
{
  connect(m_combobox, SIGNAL(activated(int)), this, SLOT(stateChanged(int)));
  connect(m_setting, SIGNAL(changed(bool)), this, SLOT(setState()));
  setState();
}

void ComboboxHandler::stateChanged(int state)
{
  m_setting->setValue(state != 0);
}

void ComboboxHandler::setState()
{
  m_combobox->setCurrentIndex(m_setting->value() ? 1 : 0);
}

SettingsDock::SettingsDock(QWidget *parent):
  QWidget(parent)
{
  m_ui.setupUi(this);

  connect(m_ui.stateSize, SIGNAL(valueChanged(int)), this, SLOT(stateSizeChanged(int)));
  connect(&Settings::instance().stateSize, SIGNAL(changed(float)), this, SLOT(setStateSize(float)));
  setStateSize(Settings::instance().stateSize.value());

  connect(m_ui.clusterHeight, SIGNAL(valueChanged(int)), this, SLOT(clusterHeightChanged(int)));
  connect(&Settings::instance().clusterHeight, SIGNAL(changed(float)), this, SLOT(setClusterHeight(float)));
  setClusterHeight(Settings::instance().clusterHeight.value());

  setupSpinbox(m_ui.branchRotation, Settings::instance().branchRotation);
  setupSpinbox(m_ui.branchTilt, Settings::instance().branchTilt);

  // new ComboboxHandler(m_ui.stateRanking, m_settings->stateRankStyleCyclic);
  // new ComboboxHandler(m_ui.clusterPositioning, m_settings->fsmStyle);
  // new ComboboxHandler(m_ui.statePositioning, m_settings->statePosStyleMultiPass);
  // new ComboboxHandler(m_ui.visualizationStyle, m_settings->clusterVisStyleTubes);
}

void SettingsDock::stateSizeChanged(int value)
{
  Settings::instance().stateSize.setValue(value / 10.0f);
}

void SettingsDock::setStateSize(float value)
{
  m_ui.stateSize->setValue((int)(value * 10.0f));
}

void SettingsDock::clusterHeightChanged(int value)
{
  Settings::instance().clusterHeight.setValue(value / 10.0f);
}

void SettingsDock::setClusterHeight(float value)
{
  m_ui.clusterHeight->setValue((int)(value * 10.0f));
}

void SettingsDock::setupSpinbox(QSpinBox *spinbox, Settings::SettingInt &setting)
{
  connect(spinbox, SIGNAL(valueChanged(int)), &setting, SLOT(setValue(int)));
  connect(&setting, SIGNAL(changed(int)), spinbox, SLOT(setValue(int)));
  spinbox->setValue(setting.value());
}
