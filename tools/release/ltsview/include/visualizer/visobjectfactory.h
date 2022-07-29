// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file visobjectfactory.h
/// \brief Header file for VisObjectFactory class

#ifndef VISOBJECTFACTORY_H
#define VISOBJECTFACTORY_H

#include <vector>
#include <QColor>
#include "vectors.h"

class PrimitiveFactory;
class VisObject;

class VisObjectFactory
{
  public:
    VisObjectFactory();
    ~VisObjectFactory();
    void clear();
    void drawObjects(PrimitiveFactory* pf,unsigned char alpha,bool texture);
    int makeObject(int primitive, std::vector<int> &ids);
    void sortObjects(const QVector3D& viewpoint);
    void updateObjectColor(int obj, QColor color);
    void updateObjectTexture(int obj, std::vector<QColor> &texColours);
    void updateObjectMatrix(int obj);
  private:
    std::vector<VisObject*> objects;
    std::vector<VisObject*> objects_sorted;
};
#endif
