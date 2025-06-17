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
        Raytracer(const Scene&);
        void createImage(const std::string&);
        Vector3d colorSet(const Ray&);
        bool isHit(const Ray&, Surface*&, Hit&);
        void worldSpace();
        bool shadowTest(Ray&, double);
        Vector3d localLight(const Ray&,  Surface*, Hit&);
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