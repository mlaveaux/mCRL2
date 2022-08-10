#ifndef MCRL2_OPENGL_CALCARC_H
#define MCRL2_OPENGL_CALCARC_H

#include <array>
#include <QVector3D>
/**
 * \brief Calculates control points for the arc described by a origin, handle
 * and destination node. \param from Position of the starting point of the arc.
 * \param via Position of the handle associated with the arc.
 * \param to Position of the end point of the arc.
 * \param selfLoop True if the arc represents a self loop.
 * \returns Four points describing a cubic Bezier curve (start, control1,
 * control2, end).
 */
static constexpr std::array<QVector3D, 4> calculateArc(const QVector3D& from,
                                                       const QVector3D& via,
                                                       const QVector3D& to,
                                                       bool selfLoop)
{
  // Pick a point a bit further from the middle point between the nodes.
  // This is an affine combination of the points 'via' and '(from + to) / 2.0f'.
  const QVector3D base = via * 1.33333f - (from + to) / 6.0f;

  if (selfLoop)
  {
    // For self-loops, the control points need to lie apart, we'll spread
    // them in x-y direction.
    const QVector3D diff =
        QVector3D::crossProduct(base - from, QVector3D(0, 0, 1));
    const QVector3D n_diff =
        diff * ((via - from).length() / (diff.length() * 2.0f));
    return std::array<QVector3D, 4>{from, base + n_diff, base - n_diff, to};
  }
  else
  {
    // Standard case: use the same position for both points.
    // return std::array<QVector3D, 4>{from, base, base, to};
    // Method 2: project the quadratic bezier curve going through the handle
    // onto a cubic bezier curve.
    const QVector3D control = via + (via - ((from + to) / 2.0f));
    return std::array<QVector3D, 4>{from, 0.33333f * from + 0.66666f * control,
                                    0.33333f * to + 0.66666f * control, to};
  }
}
#endif