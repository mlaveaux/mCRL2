// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef VECTORS_H
#define VECTORS_H

#include <QVector2D>
#include <QVector3D>

namespace Vectors {

void toPolar(float &angle, float &radius, const QVector2D& v);

QVector2D fromPolar(float angle, float radius);

}

#endif
