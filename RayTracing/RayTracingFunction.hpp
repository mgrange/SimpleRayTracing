#ifndef _RAYTRACING_H
#define _RAYTRACING_H

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
#include "RayTracingStruct.hpp"

#define NBREFLECT 10
#define NBREFRACT 10
#define ANTIALIAS 4

/*---------------------------------------------*/
/*                                             */
/*             structures                      */
/*                                             */
/*---------------------------------------------*/

class listscene{
    private:
        std::vector<Plan> allplan;
        std::vector<Sphere> allsphere;
        std::vector<Triangle> alltriangle;
        std::vector<PlanDam> alldam;
        std::vector<SphereVerre> allverre;

    public:
        listscene(){};
        ~listscene(){};

        unsigned int get_nb_plan(){
            return allplan.size();
        };
        unsigned int get_nb_sphere(){
            return allsphere.size();
        };
        unsigned int get_nb_triangle(){
            return alltriangle.size();
        };
        unsigned int get_nb_dam(){
            return alldam.size();
        };
        unsigned int get_nb_verre(){
            return allverre.size();
        };
        const Plan & get_plan(unsigned int i){
            return allplan[i];
        };
        const Sphere & get_sphere(unsigned int i){
            return allsphere[i];
        };
        const Triangle & get_triangle(unsigned int i){
            return alltriangle[i];
        };
        const PlanDam & get_dam(unsigned int i){
            return alldam[i];
        };
        const SphereVerre & get_verre(unsigned int i){
            return allverre[i];
        };
        void push(Sphere s){
            allsphere.push_back(s);
        };
        void push(Plan p){
            allplan.push_back(p);
        };
        void push(Triangle t){
            alltriangle.push_back(t);
        };
        void push(PlanDam d){
            alldam.push_back(d);
        };
        void push(SphereVerre sv){
            allverre.push_back(sv);
        };
};

class light_and_cam{
    private:
        std::vector<Point> alllight;
        Point pos_cam;
    public:
        light_and_cam(){};
        ~light_and_cam(){};
        void add_light(Point p){
            alllight.push_back(p);
        };
        Point & getlight(unsigned int i){
            return alllight[i];
        };
        Point & getcam(){
            return pos_cam;
        };
        unsigned int getnblight(){
            return alllight.size();
        };
        void setcam(Point p){
            pos_cam = p;
        };

};

/*---------------------------------------------*/
/*                                             */
/*             intersections                   */
/*                                             */
/*---------------------------------------------*/

bool intersect(const Plan& , const Ray& , Hit& , int &, listscene & allelement, light_and_cam & lumiere);
bool intersect(const PlanDam& , const Ray& , Hit& , int &, listscene & allelement, light_and_cam & lumiere);
bool intersect(const Sphere& , const Ray& , Hit& , int &, listscene & allelement, light_and_cam & lumiere);
bool intersect(const SphereVerre& , const Ray& , Hit& , int &, listscene & allelement, light_and_cam & lumiere);
bool intersect(const Triangle& , const Ray& , Hit& , int &, listscene & allelement, light_and_cam & lumiere);
bool intersect(const Ray& , Hit& , int &, listscene & allelement, light_and_cam & lumiere);

bool shadows(const Point & ,const Vector & , const Hit&,  listscene & allelement, light_and_cam & lumiere);
float speculare_light(const Vector & , const Vector & , Hit&);
float lightsource(const Point& , const Point& , Hit&, listscene & allelement, light_and_cam & lumiere);

Color antiAliasing(const unsigned int x, const unsigned int y, const Ray & ray, const Point & d0, const Vector & dx0, const Vector & dy0, const Hit & hit,  listscene & allelement, light_and_cam & lumiere);

#endif
