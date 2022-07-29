// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "settingsdialog.h"
#include <QColorDialog>

ColorButtonHandler::ColorButtonHandler(QPushButton *button, Settings::SettingColor &setting):
  QObject(button),
  m_button(button),
  m_setting(&setting)
{
  connect(m_button, SIGNAL(clicked()), this, SLOT(clicked()));
  connect(m_setting, SIGNAL(changed(QColor)), this, SLOT(setColor()));
  setColor();
}

void ColorButtonHandler::clicked()
{
  QColor color = QColorDialog::getColor(m_setting->value(), m_button->parentWidget());
  if (color.isValid())
  {
    m_setting->setValue(color);
  }
}

void ColorButtonHandler::setColor()
{
  m_button->setAutoFillBackground(true);
  QColor color = m_setting->value();
  QString colorString = QString("rgb(") + QString::number(color.red()) + ", " + QString::number(color.green()) + ", " + QString::number(color.blue()) + ")";
  m_button->setStyleSheet(QString("background-color: ") + colorString + "; color: " + colorString + ";");
}

SettingsDialog::SettingsDialog(QWidget *parent):
  QDialog(parent)
{
  m_ui.setupUi(this);

  setupSpinbox(m_ui.transparency, Settings::instance().transparency);
  new ColorButtonHandler(m_ui.backgroundColor, Settings::instance().backgroundColor);
  new ColorButtonHandler(m_ui.stateColor, Settings::instance().stateColor);
  new ColorButtonHandler(m_ui.transitionColor, Settings::instance().downEdgeColor);
  
  new ColorButtonHandler(m_ui.backpointerColorTo, Settings::instance().upEdgeColorTo);
  new ColorButtonHandler(m_ui.backpointerColorFrom, Settings::instance().upEdgeColorFrom);
  setupCheckbox(m_ui.longInterpolation_transition,
                Settings::instance().longInterpolationUpEdge);

  new ColorButtonHandler(m_ui.markColor, Settings::instance().markedColor);

  new ColorButtonHandler(m_ui.clusterColorTop, Settings::instance().clusterColorTop);
  new ColorButtonHandler(m_ui.clusterColorBottom, Settings::instance().clusterColorBottom);
  setupCheckbox(m_ui.longInterpolation, Settings::instance().longInterpolationCluster);

  new ColorButtonHandler(m_ui.simulationHistoryColor, Settings::instance().simPrevColor);
  new ColorButtonHandler(m_ui.simulationCurrentStateColor, Settings::instance().simCurrColor);
  new ColorButtonHandler(m_ui.simulationSelectedColor, Settings::instance().simSelColor);
  new ColorButtonHandler(m_ui.simulationNextStateColor, Settings::instance().simPosColor);
  setupCheckbox(m_ui.navShowBackpointers, Settings::instance().navShowBackpointers);
  setupCheckbox(m_ui.navShowStates, Settings::instance().navShowStates);
  setupCheckbox(m_ui.navShowTransitions, Settings::instance().navShowTransitions);
  setupCheckbox(m_ui.navSmoothShading, Settings::instance().navSmoothShading);
  setupCheckbox(m_ui.navLighting, Settings::instance().navLighting);
  setupCheckbox(m_ui.navTransparency, Settings::instance().navTransparency);
}

void SettingsDialog::setupSpinbox(QSpinBox *spinbox, Settings::SettingInt &setting)
{
  connect(spinbox, SIGNAL(valueChanged(int)), &setting, SLOT(setValue(int)));
  connect(&setting, SIGNAL(changed(int)), spinbox, SLOT(setValue(int)));
  spinbox->setValue(setting.value());
}

void SettingsDialog::setupCheckbox(QCheckBox *checkbox, Settings::SettingBool &setting)
{
  connect(checkbox, SIGNAL(toggled(bool)), &setting, SLOT(setValue(bool)));
  connect(&setting, SIGNAL(changed(bool)), checkbox, SLOT(setChecked(bool)));
  checkbox->setChecked(setting.value());
}
