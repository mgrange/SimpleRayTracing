#pragma once

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

bool intersect(const Plan& , const Ray& , Hit& , int &);
bool intersect(const PlanDam& , const Ray& , Hit& , int &);
bool intersect(const Sphere& , const Ray& , Hit& , int &);
bool intersect(const SphereVerre& , const Ray& , Hit& , int &);
bool intersect(const Triangle& , const Ray& , Hit& , int &);
bool intersect(const Ray& , Hit& , int &);

bool shadows(const Point & ,const Vector & , const Hit&);
float speculare_light(const Vector & , const Vector & , Hit&);
float lightsource(const Point& , const Point& , Hit&);

Color antiAliasing(const unsigned int x, const unsigned int y, const Ray & ray, const Point & d0, const Vector & dx0, const Vector & dy0, const Hit & hit);