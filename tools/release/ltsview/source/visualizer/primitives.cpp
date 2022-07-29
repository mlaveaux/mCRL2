// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file primitives.cpp
/// \brief Source file for primitive classes

#include "primitives.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include "mathutils.h"

#include "mcrl2/utilities/logger.h"

using namespace MathUtils;

/* -------- P_Sphere -------------------------------------------------------- */

P_Sphere::P_Sphere()
{
  disp_list = 0;
}

P_Sphere::~P_Sphere()
{
  glDeleteLists(disp_list,1);
}

void P_Sphere::draw()
{
  glCallList(disp_list);
}

void P_Sphere::reshape(int N,float* coss,float* sins)
{
  assert(N>0);
  int i,j,k,l;
  std::vector<GLfloat> vertices(3*((N-1)*N+2));
  std::vector<GLfloat> texCoords(((N-1)*N+2));

  vertices[0] = 0;
  vertices[1] = 0;
  vertices[2] = -1;
  texCoords[0] = 0;

  k = 3;
  l = 1;

  for (j=1; j<N; ++j)
  {
    for (i=0; i<2*N; i+=2)
    {
      vertices[k]   = sins[j] * coss[i];
      vertices[k+1] = sins[j] * sins[i];
      vertices[k+2] = -coss[j];
      texCoords[l] = sins[j] * coss[i];
      k += 3;
      ++l;
    }
  }
  vertices[k]   = 0;
  vertices[k+1] = 0;
  vertices[k+2] = 1;
  texCoords[l] = 0;

  int M = N+2;
  std::vector<GLuint> is_bot(M);

  is_bot[0] = 0;
  is_bot[1] = 1;


  for (i=N; i>=1; --i)
  {
    is_bot[M-i] = i;

  }

  
  std::vector<GLuint> is_mid((N-2)*(2*N+2));
  i = 0;
  for (j=N+1; j<=(N-2)*N+1; j+=N)
  {
    for (k=0; k<N; ++k)
    {
      is_mid[i] = j + k;
      is_mid[i+1] = is_mid[i] - N;
      i += 2;
    }
    is_mid[i] = j;
    is_mid[i+1] = is_mid[i] - N;
    i += 2;
  }

  std::vector<GLuint> is_top(M);
  j = (N-2)*N;
  is_top[0] = j + N + 1;
  i = 1;
  for (k=j+1; k<=j+N; ++k)
  {
    is_top[i] = k;

    ++i;
  }
  is_top[i] = j + 1;

  if (disp_list == 0)
  {
    disp_list = glGenLists(1);
  }
    
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);

  glVertexPointer(3,GL_FLOAT,0,vertices.data());
  glNormalPointer(GL_FLOAT,0,vertices.data());
  glTexCoordPointer(1, GL_FLOAT, 0, texCoords.data());

  glNewList(disp_list,GL_COMPILE);
  //glDrawElements(GL_TRIANGLE_FAN,M,GL_UNSIGNED_INT,is_bot.data());
  //glDrawElements(GL_TRIANGLE_FAN,M,GL_UNSIGNED_INT,is_top.data());

  M += N;
  for (int i=0; i<N-2; ++i)
  {
    //glDrawElements(GL_QUAD_STRIP,M,GL_UNSIGNED_INT,is_mid.data() + i*M);
  }
  glEndList();
}

/* -------- P_SimpleSphere -------------------------------------------------- */

P_SimpleSphere::P_SimpleSphere()
{
  GLfloat S = GLfloat(sin(PI/4));
  GLfloat C = GLfloat(cos(PI/4));
  GLfloat vertices[] = { 0,  0, -1,
                         S,  0, -C,
                         0, -S, -C,
                         -S,  0, -C,
                         0,  S,  -C,
                         1,  0,   0,
                         0, -1,   0,
                         -1,  0,   0,
                         0,  1,   0,
                         S,  0,   C,
                         0, -S,   C,
                         -S,  0,   C,
                         0,  S,   C,
                         0,  0,   1
                       };

  GLfloat texCoords[] = {0.0,
                         S,
                         0.0,
                         -S,
                         0.0,
                         1.0,
                         0.0,
                         -1.0,
                         0.0,
                         S,
                         0.0,
                         -S,
                         0.0,
                         0.0
                        };

  GLuint is_bot[] = { 0,1,2,3,4,1 };
  GLuint is_mid1[] = { 5,1,8,4,7,3,6,2,5,1 };
  GLuint is_mid2[] = { 9,5,12,8,11,7,10,6,9,5 };
  GLuint is_top[] = { 13,9,12,11,10,9 };

  disp_list = glGenLists(1);
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);

  glVertexPointer(3,GL_FLOAT,0,vertices);
  glNormalPointer(GL_FLOAT,0,vertices);
  glTexCoordPointer(1, GL_FLOAT, 0, texCoords);

  glNewList(disp_list,GL_COMPILE);
  glDrawElements(GL_TRIANGLE_FAN,6,GL_UNSIGNED_INT,is_bot);
  glDrawElements(GL_QUAD_STRIP,10,GL_UNSIGNED_INT,is_mid1);
  glDrawElements(GL_QUAD_STRIP,10,GL_UNSIGNED_INT,is_mid2);
  glDrawElements(GL_TRIANGLE_FAN,6,GL_UNSIGNED_INT,is_top);
  glEndList();
}

P_SimpleSphere::~P_SimpleSphere()
{
  glDeleteLists(disp_list,1);
}

void P_SimpleSphere::draw()
{
  glCallList(disp_list);
}

void P_SimpleSphere::reshape(int /*N*/,float* /* coss */,float* /* sins */)
{}

/* -------- P_Hemisphere ---------------------------------------------------- */

P_Hemisphere::P_Hemisphere()
{
  disp_list = 0;
}

P_Hemisphere::~P_Hemisphere()
{
  glDeleteLists(disp_list,1);
}

void P_Hemisphere::draw()
{
  glCallList(disp_list);
}

void P_Hemisphere::reshape(int N,float* coss,float* sins)
{
  assert(N>0);
  int Ndiv2 = N / 2;
  int i,j,k,l;
  std::vector<GLfloat> vertices(3 * (N * Ndiv2 + 1));
  std::vector<GLfloat> texCoords((N*Ndiv2 + 2));

  k = 0;
  l = 0;
  for (j = 0; j < Ndiv2; ++j)
  {
    for (i = 0; i < 2*N; i += 2)
    {
      vertices[k]   = coss[j] * coss[i];
      vertices[k+1] = coss[j] * sins[i];
      vertices[k+2] = sins[j];
      texCoords[l] = coss[j] * coss[i];
      k += 3;
      ++l;
    }
  }
  vertices[k]   = 0;
  vertices[k+1] = 0;
  vertices[k+2] = 1;
  texCoords[l] = 0;

  std::vector<GLuint> is_mid(N*N+N);
  i = 0;
  for (j = N; j < N*Ndiv2; j += N)
  {
    for (k = 0; k < N; ++k)
    {
      is_mid[i] = j + k;
      is_mid[i+1] = is_mid[i] - N;
      i += 2;
    }
    is_mid[i] = j;
    is_mid[i+1] = is_mid[i] - N;
    i += 2;
  }
  
  std::vector<GLuint> is_top(N+2);
  is_top[0] = N*Ndiv2;
  i = 1;
  j = (Ndiv2-1)*N;
  for (k = j; k < j+N; ++k)
  {
    is_top[i] = k;
    ++i;
  }
  is_top[i] = j;

  if (disp_list == 0)
  {
    disp_list = glGenLists(1);
  }
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glVertexPointer(3,GL_FLOAT,0,vertices.data());
  glNormalPointer(GL_FLOAT,0,vertices.data());
  glTexCoordPointer(1, GL_FLOAT, 0, texCoords.data());
  int M = 2*N + 2;
  glNewList(disp_list,GL_COMPILE);
  glDrawElements(GL_TRIANGLE_FAN,N+2,GL_UNSIGNED_INT,is_top.data());
  for (i = 0; i < Ndiv2-1; ++i)
  {
    glDrawElements(GL_QUAD_STRIP,M,GL_UNSIGNED_INT,is_mid.data() + i * M);
  }
  glEndList();
}

/* -------- P_Disc ------------------------------------------------------ */

P_Disc::P_Disc()
{
  disp_list = 0;
}

P_Disc::~P_Disc()
{
  glDeleteLists(disp_list,1);
}

void P_Disc::draw()
{
  glCallList(disp_list);
}

void P_Disc::reshape(int N,float* coss,float* sins)
{
  assert(N>0);
  GLfloat* vertices = (GLfloat*)malloc(3*N*sizeof(GLfloat));
  GLfloat* texCoords = (GLfloat*)malloc(N*sizeof(GLfloat));
  int i,j, k;
  j = 0;
  k = 0;
  for (i=0; i<2*N; i+=2)
  {
    vertices[j] = coss[i];
    vertices[j+1] = sins[i];
    vertices[j+2] = 0;
    texCoords[k] = coss[i];
    j += 3;
    ++k;
  }

  if (disp_list == 0)
  {
    disp_list = glGenLists(1);
  }

  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  glVertexPointer(3,GL_FLOAT,0,vertices);
  glTexCoordPointer(1,GL_FLOAT,0,texCoords);

  glNewList(disp_list,GL_COMPILE);
  glNormal3f(0.0f,0.0f,1.0f);
  glDrawArrays(GL_POLYGON,0,N);
  glEndList();

  free(vertices);
  free(texCoords);
}

/* -------- P_Ring ------------------------------------------------------ */

P_Ring::P_Ring(float r)
{
  disp_list = 0;
  r_top = r;
}

P_Ring::~P_Ring()
{
  glDeleteLists(disp_list,1);
}

void P_Ring::draw()
{
  glCallList(disp_list);
}

void P_Ring::reshape(int N,float* coss,float* sins)
{
  assert(N>0);
  int N3 = 3*N;
  int i,j,k,l,m;
  std::vector<GLfloat> vertices(2*N3);
  std::vector<GLfloat> texCoords(2*N+2);

  k = 0;
  l = 0;
  j = N3;
  m = N;

  for (i=0; i<2*N; i+=2)
  {
    vertices[k] = coss[i];
    vertices[k+1] = sins[i];
    vertices[k+2] = -0.5f;
    texCoords[l] = coss[i];

    k += 3;
    ++l;

    vertices[j] = r_top*coss[i];
    vertices[j+1] = r_top*sins[i];
    vertices[j+2] = 0.5f;
    texCoords[m] = r_top * coss[i];
    j += 3;
    ++m;
  }

  GLfloat* normals = (GLfloat*)malloc(2*N3*sizeof(GLfloat));
  float nz = 1.0f - r_top;
  float r = sqrt(1.0f + nz*nz);
  float nx = 1.0f / r;
  nz = nz / r;
  k = 0;
  for (i = 0; i < 2*N; i += 2)
  {
    normals[k] = coss[i] * nx;
    normals[k+1] = sins[i] * nx;
    normals[k+2] = nz;
    k += 3;
  }
  memcpy(normals+N3,normals,N3*sizeof(float));

  GLuint* is = (GLuint*)malloc((2*N+2)*sizeof(GLuint));
  k = 0;
  for (i=0; i<N; ++i)
  {
    is[k] = i + N;
    ++k;
    is[k] = i;
    ++k;

  }
  is[k] = N;
  is[k+1] = 0;

  if (disp_list == 0)
  {
    disp_list = glGenLists(1);
  }

  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);

  glVertexPointer(3,GL_FLOAT,0,vertices.data());
  glNormalPointer(GL_FLOAT,0,normals);
  glTexCoordPointer(1, GL_FLOAT, 0, texCoords.data());
  glNewList(disp_list,GL_COMPILE);
  glDrawElements(GL_QUAD_STRIP,2*N+2,GL_UNSIGNED_INT,is);
  glEndList();
  free(normals);
  free(is);
}

float P_Ring::getTopRadius()
{
  return r_top;
}

/* -------- P_TruncatedCone ------------------------------------------------------ */

P_TruncatedCone::P_TruncatedCone(P_Ring* a_ring,P_Disc* a_disc,bool t,bool b)
{
  ring = a_ring;
  disc = a_disc;
  top = t;
  bot = b;
}

P_TruncatedCone::~P_TruncatedCone()
{
}

void P_TruncatedCone::draw()
{
  if (bot)
  {
    glPushMatrix();
    glTranslatef(0.0f,0.0f,-0.5f);
    glRotatef(180,1.0f,0.0f,0.0f);
    disc->draw();
    glPopMatrix();
  }
  if (top)
  {
    glPushMatrix();
    glTranslatef(0.0f,0.0f,0.5f);
    glScalef(ring->getTopRadius(),ring->getTopRadius(),1.0f);
    disc->draw();
    glPopMatrix();
  }
  ring->draw();
}

void P_TruncatedCone::reshape(int /*N*/,float* /*coss*/,float* /*sins*/)
{
}

/* -------- P_ObliqueCone -------------------------------------------------------- */

P_ObliqueCone::P_ObliqueCone(float a,float r,float s)
{
  disp_list = 0;
  alpha = a;
  radius = r;
  sign = s;
}

P_ObliqueCone::~P_ObliqueCone()
{
  glDeleteLists(disp_list,1);
}

void P_ObliqueCone::draw()
{
  glCallList(disp_list);
}

void P_ObliqueCone::reshape(int /*N*/,float* /*coss*/,float* /*sins*/)
{
}

void P_ObliqueCone::reshape(int N,float* coss,float* sins,float obt)
{
  assert(N>0);
  float a = 0.5f*static_cast<float>(PI) - alpha - sign*obt;
  float sin_a = -sign*sin(a);
  float cos_a = cos(a);
  int i,j,k;
  GLfloat* vertices = (GLfloat*)malloc(3 * (N+1) * sizeof(GLfloat));
  GLfloat* texCoords = (GLfloat*)malloc((N + 1) * sizeof(GLfloat));

  vertices[0] = 0.0f;
  vertices[1] = 0.0f;
  vertices[2] = 0.0f;
  texCoords[0] = 0.0f;

  j = 3;
  k = 1;
  for (i = 0; i < 2*N; i += 2)
  {
    vertices[j]   = radius * coss[i] * cos_a;
    vertices[j+1] = radius * sins[i];
    vertices[j+2] = 1.0f - radius * coss[i] * sin_a;
    texCoords[k] = coss[i];
    j += 3;
    ++k;
  }

  GLfloat* normals = (GLfloat*)malloc(3 * (N+1) * sizeof(GLfloat));
  float nz = -radius;
  float r = sqrt(1.0f + nz*nz);
  float nx = 1.0f / r;
  nz = nz / r;
  normals[0] = 0.0f;
  normals[1] = 0.0f;
  normals[2] = -1.0f;
  j = 3;
  for (i = 0; i < 2*N; i += 2)
  {
    normals[j]   = coss[i] * nx * cos_a - nz * sin_a;
    normals[j+1] = sins[i] * nx;
    normals[j+2] = coss[i] * nx * sin_a + nz * cos_a;
    j += 3;
  }

  GLuint* is = (GLuint*)malloc((N+2) * sizeof(GLuint));
  is[0] = 0;
  for (i = 1; i <= N; ++i)
  {
    is[i] = N+1-i;
  }
  is[N+1] = N;

  if (disp_list == 0)
  {
    disp_list = glGenLists(1);
  }

  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  glVertexPointer(3,GL_FLOAT,0,vertices);
  glNormalPointer(GL_FLOAT,0,normals);
  glTexCoordPointer(1, GL_FLOAT, 0, texCoords);
  glNewList(disp_list,GL_COMPILE);
  glDrawElements(GL_TRIANGLE_FAN,N+2,GL_UNSIGNED_INT,is);
  glEndList();
  free(vertices);
  free(texCoords);
  free(normals);
  free(is);
}
