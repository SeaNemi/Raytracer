//William Whatley

#ifndef RAYTRACER_H
#define RAYTRACER_H

//includes the other file with the objects necessary for the scene
#include "Scene.h"

//libraries included
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <limits>

const double BIAS = 1e-4;

//used to determine which object is hit
enum class Shape{
    NOTHING,
    TRIANGLE,
    SPHERE,
    POLYGON,
};
//class Raytracer used to depict the scene
class Raytracer{
    public:
        Raytracer(Scene& scene);
        void createImage(std::string output);
        Vector3d colorSet(Ray ray);
        bool isHit(Ray ray, Surface*& currSurface, Hit &hr);
        void worldSpace();
        bool shadowTest(Ray ray, double distance);
        Vector3d localLight(Ray ray,  Surface* currSurface, Hit& hit);
        bool color;
        bool phong;
        bool stratified;
        bool reflections;
        bool shadows;
        bool verbose;
        int samples;
        double aperture;
        int maxraydepth;
    
    private:
        Scene m_scene;
        Vector3d vecW;
        Vector3d vecU;
        Vector3d vecV;
        double magD;
};

#endif