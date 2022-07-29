// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "infodock.h"
#include "cluster.h"
#include <QTableWidget>
#include <QTableWidgetItem>

static inline QTableWidgetItem *item()
{
  QTableWidgetItem *item = new QTableWidgetItem();
  item->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
  return item;
}

static void initTable(QTableWidget *table)
{
  for (int i = 0; i < table->rowCount(); i++)
  {
    for (int j = 0; j < table->columnCount(); j++)
    {
      if (!table->item(i, j))
      {
        table->setItem(i, j, item());
      }
      else
      {
        table->item(i, j)->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
      }
    }
  }
}

InfoDock::InfoDock(QWidget *parent, LtsManager *ltsManager, MarkManager *markManager):
  QWidget(parent),
  m_ltsManager(ltsManager),
  m_markManager(markManager)
{
  m_ui.setupUi(this);

  initTable(m_ui.ltsTable);
  initTable(m_ui.clusterTable);
  initTable(m_ui.stateTable);

  connect(m_ltsManager, SIGNAL(ltsChanged(LTS *)), this, SLOT(ltsChanged()));
  connect(m_ltsManager, SIGNAL(ltsZoomed(LTS *)), this, SLOT(ltsChanged()));
  connect(m_markManager, SIGNAL(statisticsChanged()), this, SLOT(markStatisticsChanged()));
  connect(m_ltsManager, SIGNAL(selectionChanged()), this, SLOT(selectionChanged()));
}

void InfoDock::ltsChanged()
{
  m_ui.ltsTable->item(0, 1)->setText(QString::number(m_ltsManager->lts()->getNumStates()));
  m_ui.ltsTable->item(1, 1)->setText(QString::number(m_ltsManager->lts()->getNumTransitions()));
  m_ui.ltsTable->item(2, 1)->setText(QString::number(m_ltsManager->lts()->getNumClusters()));
  m_ui.ltsTable->item(3, 1)->setText(QString::number(m_ltsManager->lts()->getNumRanks()));
}

void InfoDock::markStatisticsChanged()
{
  m_ui.ltsTable->item(4, 1)->setText(QString::number(m_markManager->markedStates()));
  m_ui.ltsTable->item(5, 1)->setText(QString::number(m_markManager->markedTransitions()));
}

void InfoDock::selectionChanged()
{
  LTS *lts = m_ltsManager->lts();
  int parameters = lts->getNumParameters();
  m_ui.clusterTable->setRowCount(parameters + 3);
  m_ui.stateTable->setRowCount(parameters);
  Cluster *cluster = m_ltsManager->selectedCluster();
  State *state = m_ltsManager->selectedState();

  QString statesInCluster;
  if (cluster)
  {
    statesInCluster = QString::number(cluster->getNumStates());
  }
  m_ui.clusterTable->item(0, 1)->setText(statesInCluster);

  for (int i = 0; i < parameters; i++)
  {
    QString parameterName = QString::fromStdString(lts->getParameterName(i));
    QString values;
    if (cluster)
    {
      std::set<std::string> valueSet = lts->getClusterParameterValues(cluster, i);
      bool first = true;
      values += "{ ";
      for (std::set<std::string>::iterator j = valueSet.begin(); j != valueSet.end(); j++)
      {
        if (first)
        {
          first = false;
        }
        else
        {
          values += ", ";
        }
        values += QString::fromStdString(*j);
      }
      values += " }";
    }
    QString value;
    if (state)
    {
      value = QString::fromStdString(lts->getStateParameterValueStr(state, i));
    }

    m_ui.clusterTable->setItem(i + 3, 0, item());
    m_ui.clusterTable->item(i + 3, 0)->setText(parameterName);
    m_ui.clusterTable->setItem(i + 3, 1, item());
    m_ui.clusterTable->item(i + 3, 1)->setText(values);
    m_ui.stateTable->setItem(i, 0, item());
    m_ui.stateTable->item(i, 0)->setText(parameterName);
    m_ui.stateTable->setItem(i, 1, item());
    m_ui.stateTable->item(i, 1)->setText(value);
  }
}
