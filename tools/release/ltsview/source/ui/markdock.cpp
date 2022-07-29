// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "markdock.h"
#include "markstateruledialog.h"
#include <QMessageBox>

MarkDock::MarkDock(QWidget *parent, MarkManager *markManager):
  QWidget(parent),
  m_markManager(markManager)
{
  m_ui.setupUi(this);

  m_markRuleColors = QList<QColor>()
    << QColor(228, 26, 28)
    << QColor(55, 126, 184)
    << QColor(77, 175, 74)
    << QColor(152, 78, 163)
    << QColor(255, 127, 0)
    << QColor(255, 255, 51)
    << QColor(166, 86, 40)
    << QColor(247, 129, 191)
    << QColor(153, 153, 153);
  m_markRuleNextColorIndex = 0;

  connect(m_markManager, SIGNAL(ltsChanged()), this, SLOT(loadLts()));
  connect(m_ui.noMarks, SIGNAL(clicked()), this, SLOT(markStyleClicked()));
  connect(m_ui.markDeadlocks, SIGNAL(clicked()), this, SLOT(markStyleClicked()));
  connect(m_ui.markStates, SIGNAL(clicked()), this, SLOT(markStyleClicked()));
  connect(m_ui.markTransitions, SIGNAL(clicked()), this, SLOT(markStyleClicked()));
  connect(m_markManager, SIGNAL(markStyleChanged(MarkStyle)), this, SLOT(setMarkStyle(MarkStyle)));
  connect(m_ui.clusterMatchStyle, SIGNAL(activated(int)), this, SLOT(clusterMatchStyleChanged(int)));
  connect(m_markManager, SIGNAL(clusterMatchStyleChanged(MatchStyle)), this, SLOT(setClusterMatchStyle(MatchStyle)));
  connect(m_ui.stateMatchStyle, SIGNAL(activated(int)), this, SLOT(stateMatchStyleChanged(int)));
  connect(m_markManager, SIGNAL(stateMatchStyleChanged(MatchStyle)), this, SLOT(setStateMatchStyle(MatchStyle)));
  connect(m_ui.add, SIGNAL(clicked()), this, SLOT(addMarkRule()));
  connect(m_markManager, SIGNAL(markRuleAdded(MarkRuleIndex)), this, SLOT(markRuleAdded(MarkRuleIndex)));
  connect(m_ui.markRuleList, SIGNAL(itemDoubleClicked(QListWidgetItem *)), this, SLOT(editMarkRule(QListWidgetItem *)));
  connect(m_ui.markRuleList, SIGNAL(itemChanged(QListWidgetItem *)), this, SLOT(enableMarkRule(QListWidgetItem *)));
  connect(m_markManager, SIGNAL(markRuleChanged(MarkRuleIndex)), this, SLOT(markRuleChanged(MarkRuleIndex)));
  connect(m_ui.remove, SIGNAL(clicked()), this, SLOT(removeMarkRule()));
  connect(m_markManager, SIGNAL(markRuleRemoved(MarkRuleIndex)), this, SLOT(markRuleRemoved(MarkRuleIndex)));
  connect(m_ui.markedActionList, SIGNAL(itemChanged(QListWidgetItem *)), this, SLOT(actionLabelChanged(QListWidgetItem *)));
  connect(m_markManager, SIGNAL(actionMarked(int, bool)), this, SLOT(setActionMarked(int, bool)));

  setMarkStyle(m_markManager->markStyle());
  setClusterMatchStyle(m_markManager->clusterMatchStyle());
  setStateMatchStyle(m_markManager->stateMatchStyle());
  loadLts();
}

void MarkDock::loadLts()
{
  m_actions.clear();
  m_actionNumbers.clear();
  m_actionPositions.clear();
  m_ui.markedActionList->clear();
  if (m_markManager->lts())
  {
    for (int i = 0; i < m_markManager->lts()->getNumActionLabels(); i++)
    {
      QString label = QString::fromStdString(m_markManager->lts()->getActionLabel(i));
      m_actionNumbers[label] = i;
      m_actions.push_back(label);
    }
    for (QMap<QString, int>::iterator i = m_actionNumbers.begin(); i != m_actionNumbers.end(); i++)
    {
      int index = m_ui.markedActionList->count();
      m_actionPositions[i.value()] = index;
      m_ui.markedActionList->addItem(i.key());
      m_ui.markedActionList->item(index)->setFlags(Qt::ItemIsSelectable | Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
      m_ui.markedActionList->item(index)->setCheckState(Qt::Unchecked);
    }
  }
}

void MarkDock::markStyleClicked()
{
  if (m_ui.noMarks->isChecked())
  {
    m_markManager->setMarkStyle(NO_MARKS);
  }
  else if (m_ui.markDeadlocks->isChecked())
  {
    m_markManager->setMarkStyle(MARK_DEADLOCKS);
  }
  else if (m_ui.markStates->isChecked())
  {
    m_markManager->setMarkStyle(MARK_STATES);
  }
  else if (m_ui.markTransitions->isChecked())
  {
    m_markManager->setMarkStyle(MARK_TRANSITIONS);
  }
}

void MarkDock::setMarkStyle(MarkStyle style)
{
  if (style == NO_MARKS)
  {
    m_ui.noMarks->click();
  }
  else if (style == MARK_DEADLOCKS)
  {
    m_ui.markDeadlocks->click();
  }
  else if (style == MARK_STATES)
  {
    m_ui.markStates->click();
  }
  else if (style == MARK_TRANSITIONS)
  {
    m_ui.markTransitions->click();
  }
}

void MarkDock::clusterMatchStyleChanged(int index)
{
  m_markManager->setClusterMatchStyle(index == 0 ? MATCH_ANY : MATCH_ALL);
}

void MarkDock::setClusterMatchStyle(MatchStyle style)
{
  m_ui.clusterMatchStyle->setCurrentIndex(style == MATCH_ANY ? 0 : 1);
}

void MarkDock::stateMatchStyleChanged(int index)
{
  m_markManager->setStateMatchStyle(index == 0 ? MATCH_ANY : index == 1 ? MATCH_ALL : MATCH_MULTI);
}

void MarkDock::setStateMatchStyle(MatchStyle style)
{
  m_ui.stateMatchStyle->setCurrentIndex(style == MATCH_ANY ? 0 : style == MATCH_ALL ? 1 : 2);
}

void MarkDock::addMarkRule()
{
  if (m_markManager->lts() && m_markManager->lts()->getNumParameters() > 0)
  {
    QColor color = m_markRuleColors[m_markRuleNextColorIndex];
    m_markRuleNextColorIndex = (m_markRuleNextColorIndex + 1) % m_markRuleColors.size();

    MarkStateRuleDialog dialog(this, m_markManager->lts(), color, 0, false, std::set<int>());
    if (dialog.exec() == QDialog::Accepted)
    {
      MarkRule rule;
      rule.active = true;
      rule.color = dialog.color();
      rule.parameter = dialog.parameter();
      rule.negated = dialog.negated();
      rule.values = dialog.values();
      m_markManager->addMarkRule(rule);
      m_markManager->setMarkStyle(MARK_STATES);
    }
  }
  else if(m_markManager->lts())
  {
      QMessageBox msgBox(QMessageBox::Information, "LtsView", "State marking for AUT files is not supported, since they do not contained detailed state information.",
                         QMessageBox::Ok, this, Qt::WindowCloseButtonHint);
      msgBox.exec();
  }
}

void MarkDock::markRuleAdded(MarkRuleIndex index)
{
  MarkListItem *item = new MarkListItem(markRuleDescription(index), index);
  m_markListItems[index] = item;
  m_ui.markRuleList->addItem(item);
  item->setFlags(Qt::ItemIsSelectable | Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
  item->setCheckState(Qt::Checked);
}

void MarkDock::editMarkRule(QListWidgetItem *item)
{
  MarkRuleIndex index = static_cast<MarkListItem *>(item)->index;
  MarkRule rule = m_markManager->markRule(index);
  MarkStateRuleDialog dialog(this, m_markManager->lts(), rule.color, rule.parameter, rule.negated, rule.values);
  if (dialog.exec() == QDialog::Accepted)
  {
    rule.color = dialog.color();
    rule.parameter = dialog.parameter();
    rule.negated = dialog.negated();
    rule.values = dialog.values();
    m_markManager->setMarkRule(index, rule);
    m_markManager->setMarkStyle(MARK_STATES);
  }
}

void MarkDock::enableMarkRule(QListWidgetItem *item)
{
  MarkRuleIndex index = static_cast<MarkListItem *>(item)->index;
  MarkRule rule = m_markManager->markRule(index);
  rule.active = item->checkState() == Qt::Checked;
  if (rule.active != m_markManager->markRule(index).active)
  {
    m_markManager->setMarkRule(index, rule);
    m_markManager->setMarkStyle(MARK_STATES);
  }
}

void MarkDock::markRuleChanged(MarkRuleIndex index)
{
  m_markListItems[index]->setCheckState(m_markManager->markRule(index).active ? Qt::Checked : Qt::Unchecked);
  m_markListItems[index]->setText(markRuleDescription(index));
}

void MarkDock::removeMarkRule()
{
  QList<QListWidgetItem *> selection = m_ui.markRuleList->selectedItems();
  if (!selection.isEmpty())
  {
    MarkRuleIndex index = static_cast<MarkListItem *>(selection[0])->index;
    m_markManager->removeMarkRule(index);
  }
}

void MarkDock::markRuleRemoved(MarkRuleIndex index)
{
  delete m_ui.markRuleList->takeItem(m_ui.markRuleList->row(m_markListItems[index]));
  m_markListItems.remove(index);
}

void MarkDock::actionLabelChanged(QListWidgetItem *item)
{
  bool checked = item->checkState() == Qt::Checked;
  int action = m_actionNumbers[item->text()];
  if (m_markManager->isActionMarked(action) != checked)
  {
    m_markManager->setActionMarked(action, checked);
    m_markManager->setMarkStyle(MARK_TRANSITIONS);
  }
}

void MarkDock::setActionMarked(int action, bool marked)
{
  m_ui.markedActionList->item(m_actionPositions[action])->setCheckState(marked ? Qt::Checked : Qt::Unchecked);
}

QString MarkDock::markRuleDescription(MarkRuleIndex index) const
{
  MarkRule rule = m_markManager->markRule(index);
  QString output;
  output += QString::fromStdString(m_markManager->lts()->getParameterName(rule.parameter));
  output += rule.negated ? " not in { " : " in { ";
  bool first = true;
  for (const int& i: rule.values)
  {
    if (first)
    {
      first = false;
    }
    else
    {
      output += ", ";
    }
    output += QString::fromStdString(m_markManager->lts()->getParameterDomain(rule.parameter)[i]);
  }
  output += "}";
  return output;
}
