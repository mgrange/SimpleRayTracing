#ifndef _RAYTRACING_STRUCT_H
#define _RAYTRACING_STRUCT_H

#include <cstdio>
#include <cmath>
#include <iostream>

#include "vec.h"
#include "mat.h"
#include "color.h"
#include "mesh.h"
#include "wavefront.h"
#include "image.h"
#include "orbiter.h"

/*---------------------------------------------*/
/*                                             */
/*             structures                      */
/*                                             */
/*---------------------------------------------*/

struct Ray
{
    Point origin;
    Vector direction;
};

/*un plan est defini par un point et une normale*/
struct Plan
{
    Point origin;
    Vector normal;
    Color color;
    float idref;
};

/*plan avec damier*/
struct PlanDam
{
    Point origin;
    Vector normal;
    Color color1;
    Color color2;
    float idref;
};

/*une sphere est définie par un centre et un diametre ou rayon*/
struct Sphere
{
    Point origin;
    float radius;
    Color color;
    float idref;
};

/*une sphere est définie par un centre et un diametre ou rayon*/
struct SphereVerre
{
    Point origin;
    float radius;
    Color color;
    float idref;
};

/*un triangle est définie par ces 3 points du plan*/
struct Triangle
{
    Point x1;
    Point x2;
    Point x3;
    Color color;
    float idref;
};

/*Point d'intersection avec un objet de la scene*/
struct Hit
{
    Point p;
    Vector n;
    Color color;
    float t;
};

/*---------------------------------------------*/
/*                                             */
/*          fonction D'initialisation          */
/*               Des Structures                */
/*                                             */
/*---------------------------------------------*/

Ray make_ray(const Point& o, const Point& e);

Plan make_plan(const Point& o, const Vector& v, Color color);

Plan make_plan(const Point& o, const Vector& v, Color color, float idref);

PlanDam make_planDam(const Point& o, const Vector& v, Color color1, Color color2);

PlanDam make_planDam(const Point& o, const Vector& v, Color color1, Color color2, float idref);

Sphere make_sphere(const Point& o, float r, Color color);

SphereVerre make_sphereverre(const Point& o, float r, Color color);

Sphere make_sphere(const Point& o, float r, Color color, float idref);

SphereVerre make_sphereverre(const Point& o, float r, Color color, float idref);

Triangle make_triangle(const Point& x1, const Point& x2, const Point& x3, Color c);

Triangle make_triangle(const Point& x1, const Point& x2, const Point& x3, Color c, float idref);

std::vector<Triangle> make_damier(const Point& pt,const Plan& p, int l, int c);

#endif
