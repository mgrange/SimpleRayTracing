#pragma once

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

struct Triangle
{
    Point x1;
    Point x2;
    Point x3;
    Color color;
    float idref;
};

bool intersect(const Plan& , const Ray& , Hit& , int &);
bool intersect(const PlanDam& , const Ray& , Hit& , int &);
bool intersect(const Sphere& , const Ray& , Hit& , int &);
bool intersect(const Triangle& , const Ray& , Hit& , int &);
bool intersect(const Ray& , Hit& , int &);