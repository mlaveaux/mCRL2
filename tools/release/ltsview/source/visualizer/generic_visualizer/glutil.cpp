#include "glutil.h"


// Author(s): Rimco Boudewijns and Sjoerd Cranen
/// \brief Takes 3 points and returns 4 control points to be used by an ArcShader
constexpr std::array<QVector3D, 4> GlUtil::calculateArc(const QVector3D& from, const QVector3D& via,
    const QVector3D& to, bool selfLoop){
        // Pick a point a bit further from the middle point between the nodes.
        // This is an affine combination of the points 'via' and '(from + to) / 2.0f'.
        const QVector3D base = via * 1.33333f - (from + to) / 6.0f;

        if (selfLoop)
        {
            // For self-loops, the control points need to lie apart, we'll spread
            // them in x-y direction.
            const QVector3D diff = QVector3D::crossProduct(base - from, QVector3D(0, 0, 1));
            const QVector3D n_diff = diff * ((via - from).length() / (diff.length() * 2.0f));
            return std::array<QVector3D, 4>{from, base + n_diff, base - n_diff, to};
        }
        else
        {
            // Standard case: use the same position for both points.
            //return std::array<QVector3D, 4>{from, base, base, to};
            // Method 2: project the quadratic bezier curve going through the handle onto a cubic bezier curve.
            const QVector3D control = via + (via - ((from + to) / 2.0f));
            return std::array<QVector3D, 4>{from,
            0.33333f * from + 0.66666f * control,
            0.33333f * to + 0.66666f * control,
            to};
        }
    }