#include "RayTracingStruct.hpp"

/*---------------------------------------------*/
/*                                             */
/*          fonction D'initialisation          */
/*               Des Structures                */
/*                                             */
/*---------------------------------------------*/

Ray make_ray(const Point& o, const Point& e)
{
    Ray r;
    r.origin = o;
    r.direction = make_vector(o, e);
    return r;
}

Plan make_plan(const Point& o, const Vector& v, Color color)
{
    Plan p;
    p.origin = o;
    p.normal = v;
    p.color = color;
    p.idref = 0;
    return p;
}

Plan make_plan(const Point& o, const Vector& v, Color color, float idref)
{
    Plan p;
    p.origin = o;
    p.normal = v;
    p.color = color;
    p.idref = idref;
    return p;
}

PlanDam make_planDam(const Point& o, const Vector& v, Color color1, Color color2)
{
    PlanDam p;
    p.origin = o;
    p.normal = v;
    p.color1 = color1;
    p.color2 = color2;
    p.idref = 0;
    return p;
}

PlanDam make_planDam(const Point& o, const Vector& v, Color color1, Color color2, float idref)
{
    PlanDam p;
    p.origin = o;
    p.normal = v;
    p.color1 = color1;
    p.color2 = color2;
    p.idref = idref;
    return p;
}

Sphere make_sphere(const Point& o, float r, Color color)
{
    Sphere s;
    s.origin = o;
    s.radius = r;
    s.color = color;
    s.idref = 0;
    return s;
}

SphereVerre make_sphereverre(const Point& o, float r, Color color)
{
    SphereVerre s;
    s.origin = o;
    s.radius = r;
    s.color = color;
    s.idref = 1;
    return s;
}

Sphere make_sphere(const Point& o, float r, Color color, float idref)
{
    Sphere s;
    s.origin = o;
    s.radius = r;
    s.color = color;
    s.idref = idref;
    return s;
}

SphereVerre make_sphereverre(const Point& o, float r, Color color, float idref)
{
    SphereVerre s;
    s.origin = o;
    s.radius = r;
    s.color = color;
    s.idref = idref;
    return s;
}

Triangle make_triangle(const Point& x1, const Point& x2, const Point& x3, Color c)
{
    Triangle t;
    t.x1 = x1;
    t.x2 = x2;
    t.x3 = x3;
    t.color = c;
    t.idref = 0;
    return t;
}

Triangle make_triangle(const Point& x1, const Point& x2, const Point& x3, Color c, float idref)
{
    Triangle t;
    t.x1 = x1;
    t.x2 = x2;
    t.x3 = x3;
    t.color = c;
    t.idref = idref;
    return t;
}

std::vector<Triangle> make_damier(const Point& pt,const Plan& p, int l, int c){
    std::vector<Triangle> alltriangle;

    unsigned int nbTriangles = l * c * 2;
    Vector v1 = make_vector(5.0,0.0,0.0);
    Vector v2 = make_vector(0.0,5.0,0.0);

    Point p1 = pt + v1 + p.normal;
    Point p2 = pt + v2 + p.normal;
    Point p3 = pt + v1 + v2 + p.normal;

    alltriangle.push_back(make_triangle(pt, p1, p2, make_color(1,0,0)));
    alltriangle.push_back(make_triangle(p1,p2,p3, make_color(1,0,0)));

    return alltriangle;
}
