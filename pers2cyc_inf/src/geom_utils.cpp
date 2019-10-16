#include "geom_utils.h"
#include <CGAL/Cartesian.h>
// #include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>

typedef CGAL::Cartesian<Decimal> CGALK;

void getPlaneEquation(const Vector3d& p, const Vector3d& q, 
    const Vector3d& r, Vector4d* plane_eq) {
    CGALK::Point_3 _p(p[0], p[1], p[2]);
    CGALK::Point_3 _q(q[0], q[1], q[2]);
    CGALK::Point_3 _r(r[0], r[1], r[2]);
    CGALK::Plane_3 _plane(_p, _q, _r);

    (*plane_eq)[0] = _plane.a();
    (*plane_eq)[1] = _plane.b();
    (*plane_eq)[2] = _plane.c();
    (*plane_eq)[3] = _plane.d();
}
