// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "shaders.h"

#include "mcrl2/utilities/logger.h"

namespace
{

/// \brief A simple vertex shader just transforms the 3D vertices using a worldViewProjection matrix.
const char* g_vertexShader =
  "#version 330\n"

  "uniform mat4 g_worldViewProjMatrix;\n"

  "layout(location = 0) in vec3 vertex;\n"

  "void main(void)\n"
  "{\n"
  "   gl_Position = g_worldViewProjMatrix * vec4(vertex, 1.0f);\n"
  "}";

/// \brief A simple fragment shader that uses g_color as a fill color for the polygons.
const char* g_fragmentShader =
  "#version 330\n "

  "uniform vec4 g_color = vec4(1.0f, 1.0f, 1.0f, 1.0f);\n"

  "out vec4 fragColor;\n"

  "void main(void)\n"
  "{\n"
  "   fragColor = g_color;\n"
  "}";
} // unnamed namespace

bool GlobalShader::link()
{
  // Here we compile the vertex and fragment shaders and combine the results.
  if (!addShaderFromSourceCode(QOpenGLShader::Vertex, g_vertexShader))
  {
    mCRL2log(mcrl2::log::error) << log().toStdString();
    std::abort();
  }

  if (!addShaderFromSourceCode(QOpenGLShader::Fragment, g_fragmentShader))
  {
    mCRL2log(mcrl2::log::error) << log().toStdString();
    std::abort();
  }

  if (!QOpenGLShaderProgram::link())
  {
    mCRL2log(mcrl2::log::error) << "Could not link shader program:" << log().toStdString();
    std::abort();
  }

  m_worldViewProjMatrix_location = uniformLocation("g_worldViewProjMatrix");
  if (m_worldViewProjMatrix_location == -1)
  {
    mCRL2log(mcrl2::log::warning) << "The global shader has no uniform named g_worldViewProjMatrix.\n";
  }

  m_color_location = uniformLocation("g_color");
  if (m_color_location == -1)
  {
    mCRL2log(mcrl2::log::warning) << "The global shader has no uniform named g_color.\n";
  }

  return true;
}

namespace{
  const char* g_vertexShaderInstanced = 
  "#version 330\n"
  "layout(location = 0) in vec3 vertex;\n"
  "layout(location = 2) in vec4 color;\n"
  // "layout(location = 3) in mat4 MVP;\n"
  "layout(location = 3) in vec3 offset;\n"

  "out vec4 vColor;\n"

  "uniform mat4 V;\n"
  "uniform mat4 VP;\n"
  "uniform vec3 eye;\n"
  "uniform float scale;\n"

  "void main(void)\n"
  "{\n"
  "   mat3 Vt  = transpose(mat3(V));\n"
  "   vec3 localPos = (mat3(scale) * vertex);\n"
  "   vec3 pos = (Vt*localPos + offset);\n"
  "   vec3 cam = normalize(mat3(V)*vec3(0, 0, 1));\n"
  "   vec3 p2c = normalize(eye - pos);\n"
  "   float c = dot(cam, p2c);\n"
  "   float theta = acos(c);\n"
  "   float s = sin(theta);\n"
  "   vec3 r  = normalize(cross(cam, p2c));\n"
  "   localPos = mat3(V) * localPos;\n"
  "   pos = length(cross(r, localPos)) > 0.01f ? ((1-c)*dot(localPos, r)*r + c*localPos + s*cross(r, localPos)) + offset : localPos + offset;\n"
  "   gl_Position = VP * vec4(pos, 1);\n"
  "   vColor = color;\n"
  "}";

  const char* g_fragmentShaderInstanced = 
  "#version 330\n "
  "in vec4 vColor;"
  "out vec4 fragColor;\n"
  "void main(void)\n"
  "{\n"
  "   fragColor = vColor;\n"
  "}";
}

bool NodeShaderInstanced::link(){
  if (!addShaderFromSourceCode(QOpenGLShader::Vertex, g_vertexShaderInstanced))
  {
    mCRL2log(mcrl2::log::error) << log().toStdString();
    std::abort();
  }

  if (!addShaderFromSourceCode(QOpenGLShader::Fragment, g_fragmentShaderInstanced))
  {
    mCRL2log(mcrl2::log::error) << log().toStdString();
    std::abort();
  }

  if (!QOpenGLShaderProgram::link())
  {
    mCRL2log(mcrl2::log::error) << "Could not link shader program:" << log().toStdString();
    std::abort();
  }

  u_eye_loc = uniformLocation("eye");
  if (u_eye_loc == -1)
  {
    mCRL2log(mcrl2::log::warning) << "The instanced node shader has no uniform named \"eye\".\n";
  }

  u_VP_loc = uniformLocation("VP");
  if (u_VP_loc == -1){
    mCRL2log(mcrl2::log::warning) << "The instanced node shader has no uniform named \"VP\".\n";
  }

  u_V_loc = uniformLocation("V");
  if (u_V_loc == -1){
    mCRL2log(mcrl2::log::warning) << "The instanced node shader has no uniform named \"V\".\n";
  }

  u_scale_loc = uniformLocation("scale");
  if (u_scale_loc == -1){
    mCRL2log(mcrl2::log::warning) << "The instanced node shader has no uniform named \"scale\".\n";
  }

  return true;
}

namespace{
  const char* g_vertexShaderInstancedArrow = 
  "#version 330\n"
  "layout(location = 0) in vec3 vertex;\n"
  "layout(location = 1) in vec3 offset;\n"
  "layout(location = 2) in vec4 color;\n"
  "layout(location = 3) in vec3 direction;\n"

  "out vec4 vColor;\n"

  "uniform mat4 VP;\n"
  "uniform float scale;\n"

  "void main(void)\n"
  "{\n"
  "   vec3 p = scale * vertex;\n"
  "   float theta = acos(direction.x);\n"
  "   vec3 axis = normalize(cross(vec3(1, 0, 0), direction));\n"
  "   p = (1-direction.x)*dot(p, axis)*axis + direction.x*p + sin(theta)*cross(axis, p);\n"
  "   gl_Position = VP * vec4(p+offset, 1);\n"
  "   vColor = color;\n"
  "}";

  const char* g_fragmentShaderInstancedArrow = 
  "#version 330\n "
  "in vec4 vColor;"
  "out vec4 fragColor;\n"
  "void main(void)\n"
  "{\n"
  "   fragColor = vColor;\n"
  "}";
}

bool ArrowShaderInstanced::link(){
  if (!addShaderFromSourceCode(QOpenGLShader::Vertex, g_vertexShaderInstancedArrow))
  {
    mCRL2log(mcrl2::log::error) << log().toStdString();
    std::abort();
  }

  if (!addShaderFromSourceCode(QOpenGLShader::Fragment, g_fragmentShaderInstancedArrow))
  {
    mCRL2log(mcrl2::log::error) << log().toStdString();
    std::abort();
  }

  if (!QOpenGLShaderProgram::link())
  {
    mCRL2log(mcrl2::log::error) << "Could not link shader program:" << log().toStdString();
    std::abort();
  }

  u_VP_loc = uniformLocation("VP");
  if (u_VP_loc == -1){
    mCRL2log(mcrl2::log::warning) << "The instanced arrow shader has no uniform named \"VP\".\n";
  }

  u_scale_loc = uniformLocation("scale");
  if (u_scale_loc == -1){
    mCRL2log(mcrl2::log::warning) << "The instanced arrow shader has no uniform named \"scale\".\n";
  }

  return true;
}



namespace
{
  

/// \brief A vertex shader that produces a cubic Bezier curve.
const char* g_arcVertexShader =
  "#version 330\n"

  "uniform mat4 g_viewProjMatrix;\n"
  "uniform mat4 g_viewMatrix;\n"

  "uniform vec3 g_controlPoint[4];"

  "uniform float g_density = 0.0001f;"

  "layout(location = 0) in vec3 vertex;\n"

  "out float fogAmount;\n"

  "// Calculates the position on a cubic Bezier curve with 0 <= t <= 1.\n"
  "vec3 cubicBezier(float t)\n"
  "{"
  "   return pow(1 - t, 3) * g_controlPoint[0]"
  "        + 3 * pow(1 - t, 2) * t * g_controlPoint[1]"
  "        + 3 * (1 - t) * pow(t, 2) * g_controlPoint[2]"
  "        + pow(t, 3) * g_controlPoint[3];"
  "}\n"

  "void main(void)\n"
  "{\n"
  "   // Calculate the actual position of the vertex.\n"
  "   vec4 position = vec4(cubicBezier(vertex.x), 1.0f);\n"

  "   // Apply the fog calculation to the resulting position.\n"
  "   float distance = length(g_viewMatrix * position);\n"
  "   fogAmount = (1.0f - exp(-1.0f * pow(distance * g_density, 2)));\n"

  "   // Calculate the actual vertex position in clip space.\n"
  "   gl_Position = g_viewProjMatrix * position;\n"
  "}";

/// \brief A fragment shader that uses g_color as a fill color for the polygons and applies per vertex fogging.
const char* g_arcFragmentShader =
  "#version 330\n "

  "uniform vec3 g_color = vec3(1.0f, 1.0f, 1.0f);\n"

  "in float fogAmount;\n"

  "out vec4 fragColor;\n"

  "void main(void)\n"
  "{\n"
  "   fragColor = vec4(mix(g_color, vec3(1.0f, 1.0f, 1.0f), fogAmount), 1.0);\n"
  "}";
} // unnamed namespace

bool ArcShader::link()
{
  // Here we compile the vertex and fragment shaders and combine the results.
  if (!addShaderFromSourceCode(QOpenGLShader::Vertex, g_arcVertexShader))
  {
    mCRL2log(mcrl2::log::error) << log().toStdString();
    std::abort();
  }

  if (!addShaderFromSourceCode(QOpenGLShader::Fragment, g_arcFragmentShader))
  {
    mCRL2log(mcrl2::log::error) << log().toStdString();
    std::abort();
  }

  if (!QOpenGLShaderProgram::link())
  {
    mCRL2log(mcrl2::log::error) << "Could not link shader program:" << log().toStdString();
    std::abort();
  }

  m_viewProjMatrix_location = uniformLocation("g_viewProjMatrix");
  if (m_viewProjMatrix_location == -1)
  {
    mCRL2log(mcrl2::log::warning) << "The arc shader has no uniform named g_viewProjMatrix.\n";
  }

  m_viewMatrix_location = uniformLocation("g_viewMatrix");
  if (m_viewMatrix_location == -1)
  {
    mCRL2log(mcrl2::log::warning) << "The arc shader has no uniform named g_viewMatrix.\n";
  }

  m_controlPoints_location = uniformLocation("g_controlPoint");
  if (m_controlPoints_location == -1)
  {
    mCRL2log(mcrl2::log::warning) << "The arc shader has no uniform named g_controlPoint.\n";
  }

  m_color_location = uniformLocation("g_color");
  if (m_color_location == -1)
  {
    mCRL2log(mcrl2::log::warning) << "The arc shader has no uniform named g_color.\n";
  }

  m_fogdensity_location = uniformLocation("g_density");
  if (m_fogdensity_location == -1)
  {
    mCRL2log(mcrl2::log::warning) << "The arc shader has no uniform named g_density.\n";
  }

  return true;
}

namespace
{
/// \brief A vertex shader that produces a cubic Bezier curve.
const char* g_arcVertexShaderInstanced =
  "#version 330\n"

  "uniform mat4 g_viewProjMatrix;\n"
  "uniform mat4 g_viewMatrix;\n"

  "uniform float g_density = 0.0001f;"

  "layout(location = 0) in vec3 vertex;\n"
  "layout(location = 1) in vec3 ctrl1;\n"
  "layout(location = 2) in vec3 ctrl2;\n"
  "layout(location = 3) in vec3 ctrl3;\n"
  "layout(location = 4) in vec3 ctrl4;\n"
  "layout(location = 5) in vec3 color;\n"

  "out vec4 vColor;\n"

  "out float fogAmount;\n"

  "// Calculates the position on a cubic Bezier curve with 0 <= t <= 1.\n"
  "vec3 cubicBezier(float t)\n"
  "{"
  "   return pow(1 - t, 3) * ctrl1"
  "        + 3 * pow(1 - t, 2) * t * ctrl2"
  "        + 3 * (1 - t) * pow(t, 2) * ctrl3"
  "        + pow(t, 3) * ctrl4;"
  "}\n"

  "void main(void)\n"
  "{\n"
  "   // Calculate the actual position of the vertex.\n"
  "   vec4 position = vec4(cubicBezier(vertex.x), 1.0f);\n"

  "   // Apply the fog calculation to the resulting position.\n"
  "   float distance = length(g_viewMatrix * position);\n"
  "   fogAmount = (1.0f - exp(-1.0f * pow(distance * g_density, 2)));\n"

  "   // Calculate the actual vertex position in clip space.\n"
  "   gl_Position = g_viewProjMatrix * position;\n"
  "   vColor = vec4(color, 1);\n"
  "}";

/// \brief A fragment shader that uses g_color as a fill color for the polygons and applies per vertex fogging.
const char* g_arcFragmentShaderInstanced =
  "#version 330\n "

  "in float fogAmount;\n"
  "in vec4 vColor;\n"
  "out vec4 fragColor;\n"

  "void main(void)\n"
  "{\n"
  "   fragColor = mix(vColor, vec4(1), fogAmount);\n"
  "}";
} // unnamed namespace

bool ArcShaderInstanced::link()
{
  // Here we compile the vertex and fragment shaders and combine the results.
  if (!addShaderFromSourceCode(QOpenGLShader::Vertex, g_arcVertexShaderInstanced))
  {
    mCRL2log(mcrl2::log::error) << log().toStdString();
    std::abort();
  }

  if (!addShaderFromSourceCode(QOpenGLShader::Fragment, g_arcFragmentShaderInstanced))
  {
    mCRL2log(mcrl2::log::error) << log().toStdString();
    std::abort();
  }

  if (!QOpenGLShaderProgram::link())
  {
    mCRL2log(mcrl2::log::error) << "Could not link shader program:" << log().toStdString();
    std::abort();
  }

  m_viewProjMatrix_location = uniformLocation("g_viewProjMatrix");
  if (m_viewProjMatrix_location == -1)
  {
    mCRL2log(mcrl2::log::warning) << "The arc shader has no uniform named g_viewProjMatrix.\n";
  }

  m_viewMatrix_location = uniformLocation("g_viewMatrix");
  if (m_viewMatrix_location == -1)
  {
    mCRL2log(mcrl2::log::warning) << "The arc shader has no uniform named g_viewMatrix.\n";
  }

  m_fogdensity_location = uniformLocation("g_density");
  if (m_fogdensity_location == -1)
  {
    mCRL2log(mcrl2::log::warning) << "The arc shader has no uniform named g_density.\n";
  }

  return true;
}
