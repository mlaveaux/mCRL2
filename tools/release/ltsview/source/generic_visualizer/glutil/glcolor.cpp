#include "glcolor.h"
#include <iostream>
#include <functional>
using namespace GlUtil::Color;

qreal temp = 0;
inline qreal normalize(qreal x){
  qreal y = std::modf(x - (std::trunc(x) - 1), &temp);
  y += (x==temp);
  return y;
}

inline qreal clamp01(qreal x)
{
  return std::min(1., std::max(0., x));
}

qreal* data = new qreal[6]{0, 0, 0, 0, 0, 0};
/// TODO: Dont use bit mask for clarity
/// \brief Uses \param _long to choose which value to long interpolate
inline QColor _blend(QColor& from, QColor& to, float t, std::function<void(QColor&, qreal*, qreal*, qreal*)> getVals, std::function<QColor(qreal, qreal, qreal)> getCol, int _long){
    getVals(from, &data[0], &data[1], &data[2]);
    getVals(to, &data[3], &data[4], &data[5]);
    for (int i = 0; i < 6; ++i)
      data[i] = clamp01(data[i]);
    return getCol(
        normalize(data[0] + t * (data[3] - data[0] - ((_long & 1) != 0))),
        normalize(data[1] + t * (data[4] - data[1] - ((_long & 2) != 0))),
        normalize(data[2] + t * (data[5] - data[2] - ((_long & 4) != 0))));
}


template<>
QColor GlUtil::Color::blend<ColorMode::RGB>(QColor& from, QColor& to, float t, int _long)
{
  return _blend(from, to, t,
                [](QColor& col, qreal* a, qreal* b, qreal* c){ col.getRgbF(a, b, c); },
                [](qreal a, qreal b, qreal c){     return QColor::fromRgbF(a, b, c); },
                _long);
}

template<>
QColor GlUtil::Color::blend<ColorMode::HSV>(QColor& from, QColor& to, float t, int _long){
  return _blend(from, to, t,
                [](QColor& col, qreal* a, qreal* b, qreal* c){ col.getHsvF(a, b, c); },
                [](qreal a, qreal b, qreal c){     return QColor::fromHsvF(a, b, c); },
                _long);
}

template<>
QColor GlUtil::Color::blend<ColorMode::HSL>(QColor& from, QColor& to, float t, int _long){
  return _blend(from, to, t,
                [](QColor& col, qreal* a, qreal* b, qreal* c){ col.getHslF(a, b, c); },
                [](qreal a, qreal b, qreal c){     return QColor::fromHslF(a, b, c); },
                _long);
}

template<>
QColor GlUtil::Color::lerp<InterpolateMode::SHORT, ColorMode::RGB>(QColor& from, QColor& to, float t){
    return blend<ColorMode::RGB>(from, to, t, 0);
}


template <>
QColor GlUtil::Color::lerp<InterpolateMode::SHORT, ColorMode::HSL>(QColor& from,
                                                                   QColor& to,
                                                                   float t)
{
  if (from.hslHueF() < 0)
    from.setHslF(to.hslHueF(), from.hslSaturationF(), from.lightnessF());
  return blend<ColorMode::HSL>(from, to, t, 0);
}


template<>
QColor GlUtil::Color::lerp<InterpolateMode::SHORT, ColorMode::HSV>(QColor& from, QColor& to, float t){
  if (from.hsvHueF() < 0)
    from.setHsvF(to.hsvHueF(), from.hsvSaturationF(), from.valueF());
    return blend<ColorMode::HSV>(from, to, t, 0);
}

/// \brief Pretends that r, g, b \in [0, 1] \ {0, 1}
template <>
QColor GlUtil::Color::lerp<InterpolateMode::LONG, ColorMode::RGB>(QColor& from, QColor& to, float t){
    return blend<ColorMode::RGB>(from, to, t, 7);
}

/// \brief Only long interpolate on H
template <>
QColor GlUtil::Color::lerp<InterpolateMode::LONG, ColorMode::HSL>(QColor& from, QColor& to, float t){
  if (from.hslHueF() < 0)
    from.setHslF(to.hslHueF(), from.hslSaturationF(), from.lightnessF());
    return blend<ColorMode::HSL>(from, to, t, 1);
}

/// \brief Only long interpolate on H
template <>
QColor GlUtil::Color::lerp<InterpolateMode::LONG, ColorMode::HSV>(QColor& from, QColor& to, float t){
  if (from.hsvHueF() < 0)
    from.setHsvF(to.hsvHueF(), from.hsvSaturationF(), from.valueF());
    return blend<ColorMode::HSV>(from, to, t, 1);
}
