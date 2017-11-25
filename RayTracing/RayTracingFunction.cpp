#include "RayTracingFunction.hpp"

/*---------------------------------------------*/
/*                                             */
/*          intersections                      */
/*                                             */
/*---------------------------------------------*/

/*fonction d'intersection rayon / plan*/
bool intersect(const Plan& plan, const Ray& ray, Hit& hit, int & ref, listscene & allelement, light_and_cam & lumiere)
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
            if(intersect(newRay,temp,ref,allelement,lumiere)){
                hit.color = plan.color * (1 - plan.idref) + temp.color * plan.idref;
                /*parcours de toute les lumieres de la scene et application du coefficient sur la couleur*/
                unsigned int ite;
                float coeflight = 1;
                for(ite = 0; ite < lumiere.getnblight(); ite++){
                    coeflight = coeflight * lightsource(lumiere.getlight(ite),lumiere.getcam(),temp,allelement,lumiere);
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
bool intersect(const PlanDam& plan, const Ray& ray, Hit& hit, int & ref, listscene & allelement, light_and_cam & lumiere)
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
            if(intersect(newRay,temp,ref,allelement,lumiere)){
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
                    coeflight = coeflight * lightsource(lumiere.getlight(ite),lumiere.getcam(),temp,allelement,lumiere);
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
bool intersect(const Sphere& sphere, const Ray& ray, Hit& hit, int & ref, listscene & allelement, light_and_cam & lumiere)
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
                if(intersect(newRay,temp,ref,allelement,lumiere)){
                    hit.color = sphere.color * (1 - sphere.idref) + temp.color * sphere.idref;
                    /*parcours de toute les lumieres de la scene et application du coefficient sur la couleur*/
                    unsigned int ite;
                    float coeflight = 1;
                    for(ite = 0; ite < lumiere.getnblight(); ite++){
                        coeflight = coeflight * lightsource(lumiere.getlight(ite),lumiere.getcam(),temp,allelement,lumiere);
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

bool intersect(const SphereVerre& sphere, const Ray& ray, Hit& hit, int & ref, listscene & allelement, light_and_cam & lumiere){
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
                    if(intersect(newray, temp, ref,allelement,lumiere)){                        /*parcours de toute les lumieres de la scene et application du coefficient sur la couleur*/
                        unsigned int ite;
                        float coeflight = 1;
                        for(ite = 0; ite < lumiere.getnblight(); ite++){
                            coeflight = coeflight * lightsource(lumiere.getlight(ite),lumiere.getcam(),temp,allelement,lumiere);
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
                    if(intersect(newRay,temp1,ref,allelement,lumiere)){
                        tempcol = sphere.color * (1 - sphere.idref) + temp1.color * sphere.idref;
                        /*parcours de toute les lumieres de la scene et application du coefficient sur la couleur*/
                        unsigned int ite;
                        float coeflight = 1;
                        for(ite = 0; ite < lumiere.getnblight(); ite++){
                            coeflight = coeflight * lightsource(lumiere.getlight(ite),lumiere.getcam(),temp,allelement,lumiere);
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
bool intersect(const Triangle& triangle, const Ray& ray, Hit& hit, int & ref, listscene & allelement, light_and_cam & lumiere)
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
            if(intersect(newRay,temp,ref,allelement,lumiere)){
                hit.color = triangle.color * (1 - triangle.idref) + temp.color * triangle.idref;
                /*parcours de toute les lumieres de la scene et application du coefficient sur la couleur*/
                unsigned int ite;
                float coeflight = 1;
                for(ite = 0; ite < lumiere.getnblight(); ite++){
                    coeflight = coeflight * lightsource(lumiere.getlight(ite),lumiere.getcam(),temp,allelement,lumiere);
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
bool intersect(const Ray& ray, Hit& hit, int & ref, listscene & allelement, light_and_cam & lumiere)
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
        intersect(allelement.get_sphere(i), ray, hits[i], ref, allelement,lumiere);
    }
    for(int i = 0; i < nbPlans; i++)
    {
        intersect(allelement.get_plan(i), ray, hits[i + nbSpheres], ref, allelement,lumiere);
    }
    for(int i = 0; i < nbTriangles; i++)
    {
        intersect(allelement.get_triangle(i), ray, hits[i + nbSpheres + nbPlans], ref,allelement,lumiere);
    }
    for(int i = 0; i < nbDam; i++)
    {
        intersect(allelement.get_dam(i), ray, hits[i + nbSpheres + nbPlans + nbTriangles], ref,allelement,lumiere);
    }
    for(int i = 0; i < nbVerre; i++)
    {
        intersect(allelement.get_verre(i), ray, hits[i + nbSpheres + nbPlans + nbTriangles + nbDam], ref,allelement,lumiere);
    }

    for(unsigned int i = 0; i < allelement.get_nb_plan(); i++){
        Hit temphit;
        intersect(allelement.get_plan(i), ray, temphit, ref ,allelement,lumiere);
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
bool shadows(const Point & light,const Vector & direction, const Hit& hit, listscene & allelement, light_and_cam & lumiere){
    Hit temp;
    Point machin = hit.p + direction * 0.0001f;
    Ray ray = make_ray(machin, light);

    int ref = 99;

    if(intersect(ray, temp, ref,allelement,lumiere)){
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
float lightsource(const Point& light, const Point& cam, Hit& hit, listscene & allelement, light_and_cam & lumiere)
{
    float Ambiant = 0.3;
    Vector light_ray = normalize(make_vector(hit.p, light));
    Vector cam_ray = normalize(make_vector(hit.p, cam));
    float diffuse = fmax(0.0,dot(light_ray, hit.n));
    float spec = speculare_light(light_ray,cam_ray,hit);
    bool shadow = shadows(light, light_ray, hit, allelement, lumiere);
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

Color antiAliasing(const unsigned int x, const unsigned int y, const Ray & ray, const Point & d0, const Vector & dx0, const Vector & dy0, const Hit & hit,  listscene & allelement, light_and_cam & lumiere){
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
            intersect(raya,antia,refa,allelement,lumiere);
            temp = temp + antia.color;
        }

        color = temp/(ANTIALIAS+1);
    }else{
        color = temp;
    }
    /****************************************************************/
    return color;
}
