// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file visobjectfactory.cpp
/// \brief Source file for VisObjectFactory class

#include "visobjectfactory.h"
#include <algorithm>
#include <cstdlib>
#include "primitivefactory.h"

#include <QtOpenGL>


class VisObject
{
  public:
    VisObject();
    ~VisObject();
    float* getMatrixP() const;
    QColor getColor() const; //TODO This doesn't seem to be used;
    QVector3D getCoordinates() const;
    int getPrimitive() const;
    void setColor(QColor c);
    void setTextureColours(std::vector<QColor>& colours);
    void setPrimitive(int p);
    void draw(PrimitiveFactory* pf,unsigned char alpha);
    void drawWithTexture(PrimitiveFactory* pf, unsigned char alpha);
    void addIdentifier(int id);
  private:
    float* matrix;
    QColor color;
    GLuint texName;
    int numColours;
    int primitive;
    std::vector<int> identifiers;
};

VisObject::VisObject()
{
  matrix = (float*)malloc(16*sizeof(float));
  color = QColor(150, 150, 150);
  numColours = 0;
  primitive = 0;

  // Generate texture alias for this object.
  glGenTextures(1, &texName);
}

VisObject::~VisObject()
{
  glDeleteTextures(1, &texName);
  free(matrix);
}

float* VisObject::getMatrixP() const
{
  return matrix;
}

QColor VisObject::getColor() const
{
  return color;
}

int VisObject::getPrimitive() const
{
  return primitive;
}

QVector3D VisObject::getCoordinates() const
{
  return QVector3D(matrix[12], matrix[13], matrix[14]);
}

void VisObject::setColor(QColor c)
{
  color = c;
}

void VisObject::setTextureColours(std::vector<QColor>& colours)
{
  if (colours.size() > 0)
  {
    numColours = 1;
    // numColours := smallest power of 2 s.t. colours.size <= numColours,
    // since taking an NP2 greatly reduces performance greatly.
    while (static_cast<unsigned int>(numColours) < colours.size())
    {
      numColours = numColours << 1;
    }

    GLubyte* texture = (GLubyte*)malloc(4*numColours*sizeof(GLubyte));

    for (int i = 0; i < numColours; ++i)
    {
      int j = i % colours.size();
      texture[4*i]   = colours[j].red();
      texture[4*i+1] = colours[j].green();
      texture[4*i+2] = colours[j].blue();
      texture[4*i+3] = 255; // alpha value
    }

    glBindTexture(GL_TEXTURE_1D, texName);
    // GL_TEXTURE_1D is now an alias for texName
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, numColours, 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, texture);
    free(texture);
  }
  else
  {
    numColours = 0;
  }
}

void VisObject::setPrimitive(int p)
{
  primitive = p;
}

void VisObject::draw(PrimitiveFactory* pf,unsigned char alpha)
{
  glColor4ub(color.red(), color.green(), color.blue(), alpha);
  glPushMatrix();
  glMultMatrixf(matrix);
  for (std::size_t i = 0; i < identifiers.size(); ++i)
  {
    glPushName(identifiers[i]);
  }
  pf->drawPrimitive(primitive);
  for (std::size_t i = 0; i < identifiers.size(); ++i)
  {
    glPopName();
  }
  glPopMatrix();
}

void VisObject::drawWithTexture(PrimitiveFactory* pf, unsigned char alpha)
{
  if (numColours > 0)
  {
    // Recall all settings stored in texName
    glBindTexture(GL_TEXTURE_1D, texName);
    glEnable(GL_TEXTURE_1D);
  }

  draw(pf, alpha);

  if (numColours > 0)
  {
    glDisable(GL_TEXTURE_1D);
  }
}

void VisObject::addIdentifier(int id)
{
  identifiers.push_back(id);
}
/* -------------------- Distance -------------------------------------------- */

class Distance
{
  private:
    QVector3D viewpoint;
  public:
    explicit Distance(const QVector3D& vp) : viewpoint(vp) {}
    bool operator()(const VisObject* o1,const VisObject* o2) const;
};

bool Distance::operator()(const VisObject* o1, const VisObject* o2) const
{
  QVector3D d1 = o1->getCoordinates() - viewpoint;
  QVector3D d2 = o2->getCoordinates() - viewpoint;
  return d1.lengthSquared() > d2.lengthSquared();
}

/* -------------------- VisObjectFactory ------------------------------------ */

VisObjectFactory::VisObjectFactory()
{ }

VisObjectFactory::~VisObjectFactory()
{
  clear();
}

void VisObjectFactory::sortObjects(const QVector3D& viewpoint)
{
  stable_sort(objects_sorted.begin(),objects_sorted.end(),Distance(viewpoint));
}

void VisObjectFactory::drawObjects(PrimitiveFactory* pf,unsigned char alpha,
                                   bool texture)
{
  if (texture)
  {
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

    //glEnable(GL_TEXTURE_1D);
    for (unsigned int i = 0; i < objects_sorted.size(); ++i)
    {
      objects_sorted[i]->drawWithTexture(pf,alpha);
    }
    //glDisable(GL_TEXTURE_1D);
  }
  else
  {
    for (unsigned int i = 0; i < objects_sorted.size(); ++i)
    {
      objects_sorted[i]->drawWithTexture(pf, alpha);
    }
  }
}

void VisObjectFactory::clear()
{
  for (unsigned int i = 0; i < objects.size(); ++i)
  {
    delete objects[i];
  }
  objects.clear();
  objects_sorted.clear();
}

int VisObjectFactory::makeObject(int primitive, std::vector<int> &ids)
{
  VisObject* vo = new VisObject();
  glGetFloatv(GL_MODELVIEW_MATRIX,(GLfloat*)vo->getMatrixP());
  vo->setPrimitive(primitive);

  for (std::size_t i = 0; i < ids.size(); ++i)
  {
    vo->addIdentifier(ids[i]);
  }

  objects.push_back(vo);
  objects_sorted.push_back(vo);
  return static_cast<int>(objects.size())-1;
}

void VisObjectFactory::updateObjectMatrix(int obj)
{
  glGetFloatv(GL_MODELVIEW_MATRIX,(GLfloat*)objects[obj]->getMatrixP());
}

void VisObjectFactory::updateObjectColor(int obj,QColor color)
{
  objects[obj]->setColor(color);
}

void VisObjectFactory::updateObjectTexture(int obj, std::vector<QColor> &colours)
{
  objects[obj]->setTextureColours(colours);
}
