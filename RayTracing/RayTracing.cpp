#include "RayTracingStruct.hpp"
#include "RayTracingFunction.hpp"

listscene allelement;
light_and_cam lumiere;

/*---------------------------------------------*/
/*                                             */
/*                   main                      */
/*                                             */
/*---------------------------------------------*/
int main(int agc, char **argv)
{
    srand(time(NULL));
    std::cout << "démarrage du rendu 3D :" << std::endl;

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
    Image image= create_image(1920, 1080, 3, make_color(1, 0, 0));
    /*Image en 768p*/
    // Image image= create_image(1366, 768, 3, make_color(1, 0, 0));

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
    // allelement.push(sphere);
    // allelement.push(sphere2);
    // allelement.push(sphere3);
    // allelement.push(sphere4);
    // allelement.push(plandamier);
    // allelement.push(plan);
    // allelement.push(plan2);
    // allelement.push(plan3);
    // allelement.push(plan4);
    // allelement.push(plan5);
    // allelement.push(triangle1);
    // allelement.push(triangle2);
    // allelement.push(triangle3);
    // allelement.push(triangle4);
    // allelement.push(sphere5);
    // lumiere.add_light(light_source);
    // lumiere.add_light(light_source2);

    /*scene des miroirs*/
    // allelement.push(sphere6);
    // allelement.push(plan8);
    // allelement.push(plan9);
    // allelement.push(plan6);
    // allelement.push(plan7);
    // allelement.push(plan10);
    // lumiere.add_light(light_source4);

    /*scene miroirs et refraction*/
    allelement.push(plan11);
    allelement.push(sphere7);
    allelement.push(sphere8);
    allelement.push(sphere9);
    allelement.push(sphere10);
    allelement.push(sphere11);
    lumiere.add_light(light_source3);

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

        if(intersect(ray, hit,ref,allelement,lumiere))
        {
            color = antiAliasing(x,y, ray, d0, dx0, dy0, hit, allelement, lumiere);
            // color = hit.color;

            /*parcours de toute les lumieres de la scene et application du coefficient sur la couleur*/
            unsigned int ite;
            float coeflight = 1;
            for(ite = 0; ite < lumiere.getnblight(); ite++){
                coeflight = coeflight * lightsource(lumiere.getlight(ite),lumiere.getcam(),hit,allelement,lumiere);
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
