// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file primitives.h
/// \brief Header file for primitive classes

#ifndef PRIMITIVES_H
#define PRIMITIVES_H

#include <QtOpenGL>

/* Abstract base class */

class Primitive
{
  public:
    virtual ~Primitive() {}
    virtual void draw() = 0;

    /// Reshape the primitive
    /// \param[in] N The quality of the shape.
    /// \param coss
    /// \param sins
    /// \pre N > 0
    virtual void reshape(int N,float* coss,float* sins) = 0;
};

/* Concrete derived classes */

class P_Sphere : public Primitive
{
  public:
    P_Sphere();
    ~P_Sphere();
    void draw();
    void reshape(int N,float* coss,float* sins);
  private:
    GLuint disp_list;
};

class P_SimpleSphere : public Primitive
{
  public:
    P_SimpleSphere();
    ~P_SimpleSphere();
    void draw();
    void reshape(int N,float* coss,float* sins);
  private:
    GLuint disp_list;
};

class P_Hemisphere : public Primitive
{
  public:
    P_Hemisphere();
    ~P_Hemisphere();
    void draw();
    void reshape(int N,float* coss,float* sins);
  private:
    GLuint disp_list;
};

class P_Disc : public Primitive
{
  public:
    P_Disc();
    ~P_Disc();
    void draw();
    void reshape(int N,float* coss,float* sins);
  private:
    GLuint disp_list;
};

class P_Ring : public Primitive
{
  public:
    P_Ring(float r);
    ~P_Ring();
    void draw();
    void reshape(int N,float* coss,float* sins);
    float getTopRadius();
  private:
    GLuint disp_list;
    float r_top;
};

// a truncated cone is a ring with a top surface, a bottom surface or both
class P_TruncatedCone : public Primitive
{
  public:
    P_TruncatedCone(P_Ring* a_ring,P_Disc* a_disc,bool t,bool b);
    ~P_TruncatedCone();
    void draw();
    void reshape(int N,float* coss,float* sins);
  private:
    P_Ring* ring;
    P_Disc* disc;
    bool top;
    bool bot;
};

// an oblique cone is a cone of which the base is tilted. Imagine a normal cone
// with apex on the positive z-axis and base in the (x,y)-plane. Then this cone
// becomes "oblique" if we rotate the base around the y-axis over a certain
// angle.
// Creating a P_ObliqueCone with parameters a,r,s produces a cone with the apex in
// (0,0,0), base radius r, base center at (0,0,1) and base rotated around the
// line (x=0,z=1) over the following angle: PI/2 - a - s * Branch Tilt.
// Sign s is either 1 or -1.
class P_ObliqueCone : public Primitive
{
  public:
    P_ObliqueCone(float a,float r,float s);
    ~P_ObliqueCone();
    void draw();
    void reshape(int N,float* coss,float* sins);
    void reshape(int N,float* coss,float* sins,float obt);
  private:
    GLuint disp_list;
    float alpha;
    float radius;
    float sign;
};
#endif
