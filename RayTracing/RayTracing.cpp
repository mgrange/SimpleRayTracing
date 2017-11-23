#include <cstdio>
#include <cmath>
#include <iostream>
#include <time.h>
#include <sys/time.h>

#include "vec.h"
#include "mat.h"
#include "color.h"
#include "mesh.h"
#include "wavefront.h"
#include "image.h"
#include "orbiter.h"
#include "RayTracingFunction.hpp"
#include "RayTracingStruct.hpp"

/*---------------------------------------------*/
/*                                             */
/*          intersections                      */
/*                                             */
/*---------------------------------------------*/

/*fonction d'intersection rayon / plan*/
bool intersect(const Plan& plan, const Ray& ray, Hit& hit, int & ref)
{
    /*
    plan defini par point a et normale n
    rayon defini par point o et direction d
    t = [(a-o) . n ] / (d . n)
    si d . n == 0 alors rayon // au plan => pas d'intersection
    si t negatif, pas d'intersection
    sinon, intersection
    */
    float dn = dot(ray.direction, plan.normal);
    if (dn == 0) return false;
    Vector oa = plan.origin - ray.origin;
    float t = dot(oa, plan.normal) / dn;
    if (t < 0) return false;
    hit.t = t;
    hit.p = ray.origin + hit.t*ray.direction;
    hit.n = plan.normal;
    if(plan.idref == 0){
        hit.color = plan.color;
    }else{
        // on itére sur la prochaine reflexion
        if(ref < NBREFLECT){
            ref++;
            Point newPoint = hit.p + plan.normal * 0.0001f;
            Vector newVec = ray.direction - 2.0f * dot(hit.n, ray.direction) * hit.n;
            Ray newRay;
            newRay.origin = newPoint;
            newRay.direction = newVec;
            Hit temp;
            if(intersect(newRay,temp,ref)){
                hit.color = plan.color * (1 - plan.idref) + temp.color * plan.idref;
                /*parcours de toute les lumieres de la scene et application du coefficient sur la couleur*/
                unsigned int ite;
                float coeflight = 1;
                for(ite = 0; ite < lumiere.getnblight(); ite++){
                    coeflight = coeflight * lightsource(lumiere.getlight(ite),lumiere.getcam(),temp);
                }
                hit.color = hit.color * coeflight;
            }else{
                hit.color = plan.color;
            }
        }else{
            hit.color = plan.color;
        }
    }
    return true; // t >= 0
}

/*fonction d'intersection rayon / planDam*/
bool intersect(const PlanDam& plan, const Ray& ray, Hit& hit, int & ref)
{
    float dn = dot(ray.direction, plan.normal);
    if (dn == 0) return false;
    Vector oa = plan.origin - ray.origin;
    float t = dot(oa, plan.normal) / dn;
    if (t < 0) return false;
    hit.t = t;
    hit.p = ray.origin + hit.t*ray.direction;
    hit.n = plan.normal;

    Vector u,v;
    if(plan.normal.z == 0){
        u = normalize(make_vector(plan.normal.y,-plan.normal.x,0));
        v = cross(plan.normal, u);
    }else if(plan.normal.y == 0){
        u = normalize(make_vector(plan.normal.z,-plan.normal.x,0));
        v = cross(plan.normal, u);
    }else{
        u = normalize(make_vector(plan.normal.y,-plan.normal.z,0));
        v = cross(plan.normal, u);
    }

    float xp = dot(u,make_vector(hit.p.x,hit.p.y,hit.p.z));
    float yp = dot(v,make_vector(hit.p.x,hit.p.y,hit.p.z));

    if(plan.idref == 0){
        if(xp > 0 && yp > 0){
            if(
                (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                hit.color = plan.color1;
            }else{
                hit.color = plan.color2;
            }
        }else if(xp < 0 && yp < 0){
            if(
                (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                hit.color = plan.color1;
            }else{
                hit.color = plan.color2;
            }
        }else if(xp < 0 && yp > 0){
            if(
                (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                hit.color = plan.color2;
            }else{
                hit.color = plan.color1;
            }
        }else{
            if(
                (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                hit.color = plan.color2;
            }else{
                hit.color = plan.color1;
            }
        }
    }else{
        // on itére sur la prochaine reflexion
        if(ref < NBREFLECT){
            ref++;
            Point newPoint = hit.p + plan.normal * 0.0001f;
            Vector newVec = ray.direction - 2.0f * dot(hit.n, ray.direction) * hit.n;
            Ray newRay;
            newRay.origin = newPoint;
            newRay.direction = newVec;
            Hit temp;
            if(intersect(newRay,temp,ref)){
                if(xp > 0 && yp > 0){
                    if(
                        (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                        (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                        hit.color = plan.color1;
                    }else{
                        hit.color = plan.color2 * (1 - plan.idref) + temp.color * plan.idref;
                    }
                }else if(xp < 0 && yp < 0){
                    if(
                        (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                        (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                        hit.color = plan.color1;
                    }else{
                        hit.color = plan.color2 * (1 - plan.idref) + temp.color * plan.idref;
                    }
                }else if(xp < 0 && yp > 0){
                    if(
                        (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                        (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                        hit.color = plan.color2;
                    }else{
                        hit.color = plan.color1 * (1 - plan.idref) + temp.color * plan.idref;
                    }
                }else{
                    if(
                        (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                        (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                        hit.color = plan.color2;
                    }else{
                        hit.color = plan.color1 * (1 - plan.idref) + temp.color * plan.idref;
                    }
                }
                /*parcours de toute les lumieres de la scene et application du coefficient sur la couleur*/
                unsigned int ite;
                float coeflight = 1;
                for(ite = 0; ite < lumiere.getnblight(); ite++){
                    coeflight = coeflight * lightsource(lumiere.getlight(ite),lumiere.getcam(),temp);
                }
                hit.color = hit.color * coeflight;
            }else{
                if(xp > 0 && yp > 0){
                    if(
                        (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                        (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                        hit.color = plan.color1;
                    }else{
                        hit.color = plan.color2;
                    }
                }else if(xp < 0 && yp < 0){
                    if(
                        (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                        (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                        hit.color = plan.color1;
                    }else{
                        hit.color = plan.color2;
                    }
                }else if(xp < 0 && yp > 0){
                    if(
                        (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                        (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                        hit.color = plan.color2;
                    }else{
                        hit.color = plan.color1;
                    }
                }else{
                    if(
                        (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                        (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                        hit.color = plan.color2;
                    }else{
                        hit.color = plan.color1;
                    }
                }
            }
        }else{
            if(xp > 0 && yp > 0){
                if(
                    (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                    (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                    hit.color = plan.color1;
                }else{
                    hit.color = plan.color2;
                }
            }else if(xp < 0 && yp < 0){
                if(
                    (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                    (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                    hit.color = plan.color1;
                }else{
                    hit.color = plan.color2;
                }
            }else if(xp < 0 && yp > 0){
                if(
                    (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                    (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                    hit.color = plan.color2;
                }else{
                    hit.color = plan.color1;
                }
            }else{
                if(
                    (fmod(abs(xp),2) == 0 && fmod(abs(yp),2) == 0) ||
                    (fmod(abs(xp),2) > 0 && fmod(abs(yp),2) > 0)){
                    hit.color = plan.color2;
                }else{
                    hit.color = plan.color1;
                }
            }
        }
    }
    return true; // t >= 0
}

/*fonction d'intersection rayon / sphere*/
bool intersect(const Sphere& sphere, const Ray& ray, Hit& hit, int & ref)
{
    /* cf cours :
    sphere definie par point c et rayon r
    rayon defini par point o et direction d
    (d.d)*t*t + 2d.(o-c)*t + (o-c).(o-c) - r*r = 0
    => t = (-b +/- sqrt(b*b-4*a*c))/ 2*a
    = -(2d.(o-c)) +/- sqrt((2d.(o-c))*(2d.(o-c))-4*(d.d)*(o-c).(o-c) - r*r))/ 2*(d.d)
    si discriminant < 0, racines complexes donc pas d'intersection,
    sinon, 1 ou 2 intersections.
    */
    float dd = dot(ray.direction, ray.direction);
    Vector co = ray.origin - sphere.origin;
    float coco = dot(co, co);
    float dco = dot(ray.direction, co);
    float discriminant = 4 * dco*dco - 4 * dd*(coco - sphere.radius * sphere.radius);
    bool res = (discriminant >= 0);
    if (res) {
        float b = (2 * dot(ray.direction, (ray.origin - sphere.origin)));
        float t1 = (-b + sqrt(discriminant)) / (2 * dd);
        float t2 = (-b - sqrt(discriminant)) / (2 * dd);
        hit.t = fmin(t1, t2);
        hit.p = ray.origin + hit.t * ray.direction;

        hit.n = (hit.p - sphere.origin) / sphere.radius; //division par le rad pour normaliser

        if(sphere.idref == 0){
            hit.color = sphere.color;
        }else{
            if(ref < NBREFLECT){
                ref++;
                Point newPoint = hit.p + hit.n * 0.0001f;
                Vector newVec = ray.direction - 2.0f * dot(hit.n, ray.direction) * hit.n;
                Ray newRay;
                newRay.origin = newPoint;
                newRay.direction = newVec;
                Hit temp;
                if(intersect(newRay,temp,ref)){
                    hit.color = sphere.color * (1 - sphere.idref) + temp.color * sphere.idref;
                    /*parcours de toute les lumieres de la scene et application du coefficient sur la couleur*/
                    unsigned int ite;
                    float coeflight = 1;
                    for(ite = 0; ite < lumiere.getnblight(); ite++){
                        coeflight = coeflight * lightsource(lumiere.getlight(ite),lumiere.getcam(),temp);
                    }
                    hit.color = hit.color * coeflight;
                }else{
                    hit.color = sphere.color;
                }
            }else{
                hit.color = sphere.color;
            }
        }
    }

    return (res);
}

bool intersect(const SphereVerre& sphere, const Ray& ray, Hit& hit, int & ref){
    /* cf cours :
    sphere definie par point c et rayon r
    rayon defini par point o et direction d
    (d.d)*t*t + 2d.(o-c)*t + (o-c).(o-c) - r*r = 0
    => t = (-b +/- sqrt(b*b-4*a*c))/ 2*a
    = -(2d.(o-c)) +/- sqrt((2d.(o-c))*(2d.(o-c))-4*(d.d)*(o-c).(o-c) - r*r))/ 2*(d.d)
    si discriminant < 0, racines complexes donc pas d'intersection,
    sinon, 1 ou 2 intersections.
    */
    float dd = dot(ray.direction, ray.direction);
    Vector co = ray.origin - sphere.origin;
    float coco = dot(co, co);
    float dco = dot(ray.direction, co);
    float discriminant = 4 * dco*dco - 4 * dd*(coco - sphere.radius * sphere.radius);
    bool res = (discriminant >= 0);
    if (res) {
        float b = (2 * dot(ray.direction, (ray.origin - sphere.origin)));
        float t1 = (-b + sqrt(discriminant)) / (2 * dd);
        float t2 = (-b - sqrt(discriminant)) / (2 * dd);
        hit.t = fmin(t1, t2);
        hit.p = ray.origin + hit.t * ray.direction;

        hit.n = (hit.p - sphere.origin) / sphere.radius; //division par le rad pour normaliser

        if(sphere.idref == 0){
            hit.color = sphere.color;
        }else{
            if(ref < NBREFRACT){
                ref++;
                float cosi = -dot(ray.direction,hit.n);
                float etai = 1;
                float etat = sphere.idref;
                Vector n = hit.n;
                if(cosi < 0){
                    cosi = -cosi;
                }else{
                    float temp = etai;
                    etai = etat;
                    etat = temp;
                    n = -n;
                }
                float eta = etai/etat;
                float k = 1 - eta * eta * (1 - cosi * cosi);
                if(k < 0){
                    hit.color = sphere.color;
                }else{
                    Vector newvec = eta * ray.direction + (eta * cosi - sqrtf(k)) * hit.n;
                    Ray newray;
                    Point newpoint = hit.p + newvec * 0.0001f;
                    newray.origin = newpoint;
                    newray.direction = newvec;
                    Hit temp;
                    if(intersect(newray, temp, ref)){                        /*parcours de toute les lumieres de la scene et application du coefficient sur la couleur*/
                        unsigned int ite;
                        float coeflight = 1;
                        for(ite = 0; ite < lumiere.getnblight(); ite++){
                            coeflight = coeflight * lightsource(lumiere.getlight(ite),lumiere.getcam(),temp);
                        }
                        hit.color = temp.color * coeflight;
                    }else{
                        hit.color = sphere.color;
                    }
                    Point newPoint = hit.p + hit.n * 0.0001f;
                    Vector newVec = ray.direction - 2.0f * dot(hit.n, ray.direction) * hit.n;
                    Ray newRay;
                    newRay.origin = newPoint;
                    newRay.direction = newVec;
                    Hit temp1;
                    Color tempcol;
                    if(intersect(newRay,temp1,ref)){
                        tempcol = sphere.color * (1 - sphere.idref) + temp1.color * sphere.idref;
                        /*parcours de toute les lumieres de la scene et application du coefficient sur la couleur*/
                        unsigned int ite;
                        float coeflight = 1;
                        for(ite = 0; ite < lumiere.getnblight(); ite++){
                            coeflight = coeflight * lightsource(lumiere.getlight(ite),lumiere.getcam(),temp);
                        }
                        tempcol = tempcol * coeflight;
                    }else{
                        tempcol = sphere.color;
                    }
                    hit.color = hit.color * 0.75 + tempcol * 0.25;
                }
            }else{
                hit.color = sphere.color;
            }
        }
    }

    return res;
}

/*calcule d'intersection du rayon avec un triangle*/
bool intersect(const Triangle& triangle, const Ray& ray, Hit& hit, int & ref)
{
    Vector e1 = triangle.x2 - triangle.x1;
    Vector e2 = triangle.x3 - triangle.x1;
    Vector h = cross(ray.direction, e2);
    float a = dot(e1, h);
    if (std::abs(a) < 0.00001){ return false; }
    float f = 1.f / a;
    Vector s = ray.origin - triangle.x1;
    float u = f * dot(s, h);
    if (u < 0 || u > 1) { return false; }
    Vector q = cross(s, e1);
    float v = f * dot(ray.direction, q);
    if (v < 0 || (u + v) > 1){ return false; }
    float t = f * dot(e2, q);
    if (t < 0){ return false; }

    hit.t = t;
    hit.p = ray.origin + t * ray.direction;
    hit.n = normalize(cross((triangle.x3 - triangle.x1), (triangle.x2 - triangle.x1)));
    if(triangle.idref == 0){
            hit.color = triangle.color;
        }else{
        if(ref < NBREFLECT){
            ref++;
            Point newPoint = hit.p + hit.n * 0.0001f;
            Vector newVec = ray.direction - 2.0f * dot(hit.n, ray.direction) * hit.n;
            Ray newRay;
            newRay.origin = newPoint;
            newRay.direction = newVec;
            Hit temp;
            if(intersect(newRay,temp,ref)){
                hit.color = triangle.color * (1 - triangle.idref) + temp.color * triangle.idref;
                /*parcours de toute les lumieres de la scene et application du coefficient sur la couleur*/
                unsigned int ite;
                float coeflight = 1;
                for(ite = 0; ite < lumiere.getnblight(); ite++){
                    coeflight = coeflight * lightsource(lumiere.getlight(ite),lumiere.getcam(),temp);
                }
                hit.color = hit.color * coeflight;
            }else{
                hit.color = triangle.color;
            }
        }else{
            hit.color = triangle.color;
        }
    }
    return true;
}

/*
renvoie vrai si le rayon a touché un des objets
et que le point d'intersection est bien devant la camera / l'origine du rayon,
ainsi que la distance de l'objet le plus proche
*/
bool intersect(const Ray& ray, Hit& hit, int & ref)
{
    bool res = false;
    Hit hitmin;
    float t_init = 9999;

    int nbSpheres = allelement.get_nb_sphere();
    int nbPlans = allelement.get_nb_plan();
    int nbTriangles = allelement.get_nb_triangle();
    int nbDam = allelement.get_nb_dam();
    int nbVerre = allelement.get_nb_verre();

    Hit hit_tmp;
    hit_tmp.t = t_init;
    std::vector<Hit> hits;
    for(int i = 0; i< nbSpheres + nbPlans + nbTriangles + nbDam + nbVerre; i++){
        hits.push_back(hit_tmp);
    }

    for(int i = 0; i< nbSpheres; i++)
    {
        intersect(allelement.get_sphere(i), ray, hits[i], ref);
    }
    for(int i = 0; i < nbPlans; i++)
    {
        intersect(allelement.get_plan(i), ray, hits[i + nbSpheres], ref);
    }
    for(int i = 0; i < nbTriangles; i++)
    {
        intersect(allelement.get_triangle(i), ray, hits[i + nbSpheres + nbPlans], ref);
    }
    for(int i = 0; i < nbDam; i++)
    {
        intersect(allelement.get_dam(i), ray, hits[i + nbSpheres + nbPlans + nbTriangles], ref);
    }
    for(int i = 0; i < nbVerre; i++)
    {
        intersect(allelement.get_verre(i), ray, hits[i + nbSpheres + nbPlans + nbTriangles + nbDam], ref);
    }

    for(unsigned int i = 0; i < allelement.get_nb_plan(); i++){
        Hit temphit;
        intersect(allelement.get_plan(i), ray, temphit, ref);
        if(temphit.t > 0 && temphit.t < hitmin.t){
            hitmin = temphit;
            // printf("res ? = %d\n", res);
        }
    }

    Hit hit_min;
    hit_min.t = t_init;
    for(unsigned int i = 0; i < hits.size(); i++)
    {
        if(hits[i].t > 0 && hits[i].t < hit_min.t){
            hit_min = hits[i];
        }
    }
    if(hit_min.t==t_init) return false; //aucun objet touche
    hit = hit_min;
    return true;
}


/*---------------------------------------------*/
/*                                             */
/*               Lumières                      */
/*                                             */
/*---------------------------------------------*/

/*determine si l'objet est à l'ombre ou non*/
bool shadows(const Point & light,const Vector & direction, const Hit& hit){
    Hit temp;
    Point machin = hit.p + direction * 0.0001f;
    Ray ray = make_ray(machin, light);

    int ref = 99;

    if(intersect(ray, temp, ref)){
        float DistPI = pow(temp.p.x - hit.p.x,2) + pow(temp.p.y - hit.p.y,2) + pow(temp.p.z - hit.p.z,2);
        float DistLI = pow(light.x - hit.p.x,2) + pow(light.y - hit.p.y,2) + pow(light.z - hit.p.z,2);
        if(DistPI < DistLI){
            return true;
        }else{
            return false;
        }
    }else{
        return false;
    }
}

/*calcule la lumière spéculaire sur la scene*/
float speculare_light(const Vector & light_ray, const Vector & cam_ray, Hit& hit){
    if(dot(hit.n, light_ray) > 0){
        Vector H = normalize(light_ray + cam_ray);
        return pow(dot(hit.n, H),30);
    }
    else{
        return 0;
    }
}

/*choisi la couleur en fonction de l'objet toucher à modifié*/
float lightsource(const Point& light, const Point& cam, Hit& hit)
{
    float Ambiant = 0.3;
    Vector light_ray = normalize(make_vector(hit.p, light));
    Vector cam_ray = normalize(make_vector(hit.p, cam));
    float diffuse = fmax(0.0,dot(light_ray, hit.n));
    float spec = speculare_light(light_ray,cam_ray,hit);
    bool shadow = shadows(light, light_ray, hit);
    float pos_angle;
    if(!shadow){
        pos_angle = Ambiant + diffuse + spec;
    }else{
        pos_angle = Ambiant;
    }
    return pos_angle;
}

/*---------------------------------------------*/
/*                                             */
/*                 antialiasing                */
/*                                             */
/*---------------------------------------------*/

Color antiAliasing(const unsigned int x, const unsigned int y, const Ray & ray, const Point & d0, const Vector & dx0, const Vector & dy0, const Hit & hit){
    Color temp = hit.color/(ANTIALIAS+1);
    Color color;
    /*mise en place de l'anti-aliasing*/
    /****************************************************************/
    if(ANTIALIAS != 0){
        for(unsigned int z = 0; z < ANTIALIAS; z++){
            int refa = 0;
            Hit antia;
            Ray raya = ray;
            float newx = x + (float)(rand()/RAND_MAX - 0.5);
            float newy = y + (float)(rand()/RAND_MAX- 0.5);
            raya.origin = d0 + newx*dx0 + newy*dy0;
            intersect(raya,antia,refa);
            temp = temp + antia.color;
        }

        color = temp/(ANTIALIAS+1);
    }else{
        color = temp;
    }
    /****************************************************************/
    return color;
}

/*---------------------------------------------*/
/*                                             */
/*                   main                      */
/*                                             */
/*---------------------------------------------*/
int main(int agc, char **argv)
{
    srand(time(NULL));
    std::cout << "démarage du rendu 3D :" << std::endl;

    struct timeval begin, end;
    gettimeofday(&begin, NULL);
    // creer la camera
    const Point& cam_centre = make_point(0.0, 0.0, 0.0);
    const Point& light_source = make_point(5.0, 5.0, 10.0);
    const Point& light_source2 = make_point(-5.0, 5.0, 10.0);
    const Point& light_source3 = make_point(0,0,30);
    const Point& light_source4 = make_point(-0.5,-0.5,10);
    const float size = 10;
    Orbiter camera = make_orbiter_lookat(cam_centre, size);
    //Orbiter camera = make_orbiter(); //regarde 0,0,0 depuis une distance 5

    // creer l'image resultat
    /*Image en 4K*/
    // Image image= create_image(3840, 2160, 3, make_color(1, 0, 0));
    /*Image en 1080p*/
    // Image image= create_image(1920, 1080, 3, make_color(1, 0, 0));
    /*Image en 768p*/
    Image image= create_image(1366, 768, 3, make_color(1, 0, 0));

    // creer les objets de la scene
    Point o1 = make_point(0, -1, 0);
    Vector n1 = make_vector(0, 1, 0);
    PlanDam plandamier = make_planDam(o1, n1, make_color(1.0,0.2,0.2), make_color(0.2,1.0,0.2),0.1);

    Point o2 = make_point(1, 0, 0);
    float r1 = 1;
    Sphere sphere = make_sphere(o2, r1, make_color(0.2,1.0,0.2), 0.2);

    Point o3 = make_point(-2, 0, 0);
    float r2 = 1;
    Sphere sphere2 = make_sphere(o3, r2, make_color(0.2,0.2,1.0), 0.4);

    Point o4 = make_point(0, -1, -20);
    Vector n2 = make_vector(0, 0, 1);
    Plan plan = make_plan(o4, n2, make_color(0.1,0.4,0.6));

    Point o5 = make_point(5, 0, -5);
    float r3 = 1;
    Sphere sphere3 = make_sphere(o5, r3, make_color(1.0,0.2,0.2), 0.6);

    Point o6 = make_point(2.0,2.0,0.0);
    Point o7 = make_point(2.0,3.0,-1.0);
    Point o8 = make_point(2.0,1.0,-1.0);
    Point o9 = make_point(4.0,2.0,0.0);

    Triangle triangle1 = make_triangle(o6,o9,o7,make_color(0.9,0.2,0.9));
    Triangle triangle2 = make_triangle(o6,o8,o7,make_color(0.9,0.2,0.9));
    Triangle triangle3 = make_triangle(o6,o9,o8,make_color(0.9,0.2,0.9));
    Triangle triangle4 = make_triangle(o8,o9,o7,make_color(0.9,0.2,0.9));

    Point o10 = make_point(0, 1, 10);
    float r4 = 2;
    Sphere sphere4 = make_sphere(o10, r4, make_color(1.0,0.5,0.2));

    Point o11 = make_point(0, -1, 11);
    Vector n3 = make_vector(0, 0, 1);
    Plan plan2 = make_plan(o11, n3, make_color(0.9,0.2,0.2));

    Point o12 = make_point(6, 0, 0);
    Vector n4 = make_vector(-1, 0, 0);
    Plan plan3 = make_plan(o12, n4, make_color(0.2,0.2,0.9),1);

    Point o13 = make_point(-6, 0, 0);
    Vector n5 = make_vector(1, 0, 0);
    Plan plan4 = make_plan(o13, n5, make_color(0.2,0.9,0.2));

    Point o14 = make_point(0, 10, 0);
    Vector n6 = make_vector(0, -1, 0);
    Plan plan5 = make_plan(o14, n6, make_color(0.9,0.9,0.2));

    Point o15 = make_point(-3, 3, 0);
    float r5 = 1;
    SphereVerre sphere5 = make_sphereverre(o15, r5, make_color(1.0,0.2,0.2),0.2);

    Point o16 = make_point(0,0,-15);
    float r6 = 2;
    Sphere sphere6 = make_sphere(o16,r6,make_color(0.9,0.2,0.2));

    Point o17 = make_point(0,5,0);
    Vector n7 = make_vector(0,-1,0);
    Plan plan6 = make_plan(o17,n7,make_color(0.9,0.1,0.9));

    Point o18 = make_point(0,-5,0);
    Vector n8 = make_vector(0,1,0);
    Plan plan7 = make_plan(o18,n8,make_color(0.9,0.1,0.1));

    Point o19 = make_point(3,0,0);
    Vector n9 = make_vector(-1,0,0);
    Plan plan8 = make_plan(o19,n9,make_color(0.1,0.9,0.1),1);

    Point o20 = make_point(-3,0,0);
    Vector n10 = make_vector(1,0,0);
    Plan plan9 = make_plan(o20,n10,make_color(0.1,0.1,0.9),1);

    Point o21 = make_point(0, 0, -20);
    Vector n11 = make_vector(0, 0, 1);
    Plan plan10 = make_plan(o21, n11, make_color(0.9,0.9,0.1),0.2);

    PlanDam plan11 = make_planDam(o21, n11, make_color(0.9,0.9,0.2), make_color(0.2,0.2,0.2));

    Point o22 = make_point(0,0,-3);
    float r7 = 1;
    SphereVerre sphere7 = make_sphereverre(o22,r7,make_color(0.9,0.2,0.2),0.5);

    Point o23 = make_point(2,2,-3);
    float r8 = 1;
    Sphere sphere8 = make_sphere(o23,r8,make_color(0.9,0.2,0.2),0.5);

    Point o24 = make_point(-2,-2,-3);
    float r9 = 1;
    Sphere sphere9 = make_sphere(o24,r9,make_color(0.9,0.2,0.9),0.2);

    Point o25 = make_point(2,-2,-3);
    float r10 = 1;
    Sphere sphere10 = make_sphere(o25,r10,make_color(0.2,0.2,0.9),0.5);

    Point o26 = make_point(-2,2,-3);
    float r11 = 1;
    Sphere sphere11 = make_sphere(o26,r11,make_color(0.2,0.9,0.2),0.5);

    /*scene montrant tous calcule en une fois*/
    allelement.push(sphere);
    allelement.push(sphere2);
    allelement.push(sphere3);
    allelement.push(sphere4);
    allelement.push(plandamier);
    allelement.push(plan);
    allelement.push(plan2);
    allelement.push(plan3);
    allelement.push(plan4);
    allelement.push(plan5);
    allelement.push(triangle1);
    allelement.push(triangle2);
    allelement.push(triangle3);
    allelement.push(triangle4);
    allelement.push(sphere5);
    lumiere.add_light(light_source);
    lumiere.add_light(light_source2);

    /*scene des miroirs*/
    // allelement.push(sphere6);
    // allelement.push(plan8);
    // allelement.push(plan9);
    // allelement.push(plan6);
    // allelement.push(plan7);
    // allelement.push(plan10);
    // lumiere.add_light(light_source4);

    /*scene miroirs et refraction*/
    // allelement.push(plan11);
    // allelement.push(sphere7);
    // allelement.push(sphere8);
    // allelement.push(sphere9);
    // allelement.push(sphere10);
    // allelement.push(sphere11);
    // lumiere.add_light(light_source3);

    //boucler sur les pixels
    Point d0;
    Vector dx0, dy0;
    float fov = 45; //???... en degrés !
    Point o = orbiter_position(camera);
    orbiter_image_frame(camera, image.width, image.height, 0/*???*/, fov, d0, dx0, dy0);

    lumiere.setcam(o);

    #pragma omp parallel for
    for(int y= 0; y < image.height; y++)
    for(int x= 0; x < image.width; x++)
    {
        // generer l'origine et l'extremite du rayon
        //cf Orbiter.h

        Point e = d0 + x*dx0 + y*dy0;

        Ray ray= make_ray(o, e);

        // calculer les intersections
        Hit hit;
        Color color;
        int ref = 0;

        if(intersect(ray, hit,ref))
        {
            color = antiAliasing(x,y, ray, d0, dx0, dy0, hit);
            // color = hit.color;

            /*parcours de toute les lumieres de la scene et application du coefficient sur la couleur*/
            unsigned int ite;
            float coeflight = 1;
            for(ite = 0; ite < lumiere.getnblight(); ite++){
                coeflight = coeflight * lightsource(lumiere.getlight(ite),lumiere.getcam(),hit);
            }
            color = color * coeflight;
        }
        else
        {
            color = make_color(0, 0, 0);
        }
        image_set_pixel(image, x, y, color);
    }

    // enregistrer l'image
    write_image(image, "render.png");
    // nettoyage
    release_image(image);
    gettimeofday(&end, NULL);
    printf("Durée : %.3f secondes.\n", (double)((end.tv_sec  - begin.tv_sec) * 1000000u + end.tv_usec - begin.tv_usec) / 1.e6);
    return 0;

}
