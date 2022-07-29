// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file primitivefactory.h
/// \brief Header file for PrimitiveFactory class

#ifndef PRIMITIVEFACTORY_H
#define PRIMITIVEFACTORY_H

#include <vector>
#include "settings.h"

class ConeDB;
class Primitive;
class P_ObliqueCone;

class PrimitiveFactory: public QObject
{
  Q_OBJECT

  public:
    PrimitiveFactory();
    ~PrimitiveFactory();
    void drawPrimitive(int p);
    void drawSimpleSphere();
    int makeHemisphere();
    int makeObliqueCone(float a,float r,float s);
    int makeSphere();
    int makeTruncatedCone(float r,bool topClosed,bool bottomClosed);
    //int  makeTube();

  private slots:
    void qualityChanged();
    void branchTiltChanged();

  private:
    std::vector<Primitive*> primitives;
    std::vector<P_ObliqueCone*> oblq_cones;
    ConeDB*  coneDB;
    int     disc;
    int     simple_sphere;
    int     sphere;
    int     hemisphere;
    float*   cos_theta;
    float*   sin_theta;

    int  make_ring(float r);
    void make_disc();
    void make_simple_sphere();
    void update_geom_tables();
    void update_oblique_cones();
    void update_primitives();
};
#endif
