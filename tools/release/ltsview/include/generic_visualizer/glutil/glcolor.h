#ifndef MCRL2_OPENGL_COLOR_H
#define MCRL2_OPENGL_COLOR_H

// Assumed usage of Qt for colors
#include <QColor>
#include <cmath>

namespace GlUtil
{
namespace Color
{

/// TODO: Add CMYK?
enum class ColorMode{
    RGB,
    HSV,
    HSL
};

enum class InterpolateMode{
    SHORT,
    LONG
};

template <InterpolateMode imode, ColorMode cmode>
inline QColor lerp(QColor& from, QColor& to, float t);


template <ColorMode mode>
inline QColor blend(QColor& from, QColor& to, float t, int _long);


} // namespace Color
} // namespace GlUtil
#endif