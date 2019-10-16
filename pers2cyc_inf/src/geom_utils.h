#ifndef _GEOM_UTILS_H_
#define _GEOM_UTILS_H_

#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "utils.h"

void getPlaneEquation(const Vector3d& p, const Vector3d& q, 
    const Vector3d& r, Vector4d* plane_eq);

inline Decimal evalPlaneEquation(const Vector4d& plane_eq, const Vector3d& p) {
    return plane_eq[0]*p[0] + plane_eq[1]*p[1] + plane_eq[2]*p[2] + plane_eq[3];
}

// a-b in terms of angle radian, i.e., the radian to 
// sweep for when going from b to a
inline double radianSub(double a, double b) {
    if (a >= b) {
        return a-b;
    } else {
        return 2*M_PI-b+a;
    }
}

inline double radian2Degree(double r) {
    return r * 180.0 / M_PI;
}

#endif
