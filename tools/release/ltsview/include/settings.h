// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file settings.h
/// \brief Header file for Settings class

#ifndef SETTINGS_H
#define SETTINGS_H

#include "mcrl2/gui/setting.h"

class Settings
{
public:
  typedef mcrl2::gui::qt::SettingInt SettingInt;
  typedef mcrl2::gui::qt::SettingFloat SettingFloat;
  typedef mcrl2::gui::qt::SettingBool SettingBool;
  typedef mcrl2::gui::qt::SettingColor SettingColor;

  SettingFloat stateSize;
  SettingFloat clusterHeight;
  SettingInt branchRotation;
  SettingInt branchTilt;
  SettingInt quality;
  SettingInt transparency;
  SettingColor backgroundColor;
  SettingColor stateColor;
  SettingColor downEdgeColor;

  SettingColor upEdgeColorFrom;
  SettingColor upEdgeColorTo;
  SettingBool longInterpolationUpEdge;

  SettingColor markedColor;
  
  SettingColor clusterColorTop;
  SettingColor clusterColorBottom;
  SettingBool longInterpolationCluster;
  
  SettingColor simPrevColor;
  SettingColor simCurrColor;
  SettingColor simSelColor;
  SettingColor simPosColor;
  SettingBool stateRankStyleCyclic;
  SettingBool fsmStyle;
  SettingBool statePosStyleMultiPass;
  SettingBool clusterVisStyleTubes;
  SettingBool navShowBackpointers;
  SettingBool navShowStates;
  SettingBool navShowTransitions;
  SettingBool navSmoothShading;
  SettingBool navLighting;
  SettingBool navTransparency;
  SettingBool displayStates;
  SettingBool displayTransitions;
  SettingBool displayBackpointers;
  SettingBool displayWireframe;
  
private:
  Settings():
    stateSize(0.1f),
    clusterHeight(0.0f),
    branchRotation(0),
    branchTilt(30),
    quality(32),
    transparency(30),
    backgroundColor(QColor(150, 150, 150)),
    stateColor(QColor(255, 255, 255)),
    downEdgeColor(QColor(255, 255, 255)),
    upEdgeColorFrom(QColor(0, 0, 0)),
    upEdgeColorTo(QColor(0, 0, 255)),
    longInterpolationUpEdge(false),
    markedColor(QColor(255, 0, 0)),
    clusterColorTop(QColor(255, 255, 255)),
    clusterColorBottom(QColor(0, 0, 255)),
    longInterpolationCluster(false),
    simPrevColor(QColor(255, 0, 255)),
    simCurrColor(QColor(0, 255, 0)),
    simSelColor(QColor(255, 255, 0)),
    simPosColor(QColor(0, 255, 255)),
    stateRankStyleCyclic(false),
    fsmStyle(false),
    statePosStyleMultiPass(false),
    clusterVisStyleTubes(false),
    navShowBackpointers(true),
    navShowStates(true),
    navShowTransitions(true),
    navSmoothShading(true),
    navLighting(true),
    navTransparency(true),
    displayStates(false),
    displayTransitions(false),
    displayBackpointers(false),
    displayWireframe(false)
  {}

public:
  // singleton pattern
  static Settings& instance(){
    static Settings INSTANCE;
    return INSTANCE;
  }
};

#endif
