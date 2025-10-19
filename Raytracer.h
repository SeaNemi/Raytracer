//William Whatley
//Raytracer header file
#ifndef RAYTRACER_H
#define RAYTRACER_H

//includes the other file with the objects necessary for the scene
#include "Scene.h"

//libraries includedand constant
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
const double BIAS = 1e-6;

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
        bool color;
        bool phong;
        bool stratified;
        bool reflections;
        bool shadows;
        int samples;
        double aperture;
        int maxraydepth;
    
    private:
        //helper functions and member variables private
        bool isHit(Node*, const Ray&, Surface*&, Hit&, double, double&);
        void worldSpace();
        bool shadowTest(const Ray&, double);
        Vector3d localLight(const Ray&, Surface*, Hit&);
        Vector3d colorSet(const Ray&);
        Scene m_scene;
        Vector3d vecW;
        Vector3d vecU;
        Vector3d vecV;
        double magD;
};

#endif