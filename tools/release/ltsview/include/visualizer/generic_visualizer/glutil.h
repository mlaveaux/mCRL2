// Author(s): Ruben Vink.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_OPENGL_GLUTIL_H
#define MCRL2_OPENGL_GLUTIL_H

#include <array>
#include <QVector3D>

class GlUtil{
public:
  /**
  * \brief Calculates control points for the arc described by a origin, handle and destination node.
  * \param from Position of the starting point of the arc.
  * \param via Position of the handle associated with the arc.
  * \param to Position of the end point of the arc.
  * \param selfLoop True if the arc represents a self loop.
  * \returns Four points describing a cubic Bezier curve (start, control1, control2, end).
  */
  static constexpr std::array<QVector3D, 4> calculateArc(const QVector3D& from, const QVector3D& via,
    const QVector3D& to, bool selfLoop);

private:
    // make constructor private so we don't make objects of this type
    GlUtil() {}
};


#endif