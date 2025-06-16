//William Whatley
//Objects Header File

//Header guards
#ifndef OBJECTS_H
#define OBJECTS_H

//What's included
#include <iostream>
#include "matrix.h"
#include <string>
#include <vector>

class Ray{
    public:
        Ray(Vector3d e, Vector3d d, int dep = 0);
        Vector3d eye;
        Vector3d dir;
        int depth;
};

class Hit{
    public:
        Hit(){}
        Hit(double t,double a,double b,double g, Vector3d inte, Vector3d nor, Vector3d vie, int depth);
        double t, alpha, beta, gamma;
        Vector3d inter; //determines intersection
        Vector3d norm; //normalized perpendicular to intersection pt
        Vector3d view; //determines the viewing direction
        int raydepth; //determines the depth
        void transfer(double fil[8]);
        double fill[8];
};

class Surface{
    public:
        virtual ~Surface();
        virtual bool intersect(Ray ray, double &min, double &max, Hit &hr) = 0;
        virtual Surface* clone() = 0;
        double m_fill[8];
        bool isPatch;
};
//Polygon class
class Polygon: public Surface{
    public:
        Polygon(){}
        Polygon(int vertices, double fill[8], bool tri = false, bool pat = false);
        Polygon(const Polygon& poly);
        virtual ~Polygon();
        void clear();
        int getSize();
        void addVertex(double x, double y, double z);
        void pushBackVector(Vector3d*);
        Polygon& operator=(const Polygon& other);
        bool intersect(Ray ray, double &min, double &max, Hit &hr) override;
        int m_verticies;
        std::vector<Vector3d*> m_vertex;

        Surface* clone();
        bool isTriangle;
};

//patch class
//child of Polygon
class Patch: public Polygon{
    public:
        //functions
        Patch(int veritices, double fill[8], bool isTri = false, bool pat = true);
        Patch(const Patch& rhs);
        ~Patch();
        void addNormal(double x, double y, double z);
        void pushBackNorm(Vector3d* a);
        //member variables
        std::vector<Vector3d*> m_normal;
        Vector3d interpolate(Hit hit);
        Surface* clone() override;
};

//Sphere class
class Sphere: public Surface{
    public:
        //public
        Sphere(){}
        Sphere(const Sphere& rhs);
        Sphere(double x, double y, double z, double r, double fill[8], bool pat = false);
        Sphere& operator=(const Sphere& rhs);
        bool intersect(Ray ray, double &min, double &max, Hit &hr) override;
        Surface* clone();

        //member variables
        Vector3d coords;
        double radius;
};


//Light class
class Light{
    public:
        Light();
        Light(double x, double y, double z);
        //Vector3d intensity();
        Vector3d coords;

};

//Viewpoint class
class Viewpoint{
    public:
        //viewpoint constructor has nothing
        Viewpoint(){}
        Viewpoint(const Viewpoint& rhs);
        Viewpoint& operator=(const Viewpoint& rhs);

        //all setters
        void setAt(double x, double y, double z);
        void setFrom(double x, double y, double z);
        void setUp(double x, double y, double z);
        void setAngle(double a);
        void setHither(double g);
        void setRes(int x, int y);


        //doubles for the viewpoint
        Vector3d at;
        Vector3d from;
        Vector3d up;
        double angle;
        double hither;
        double res[2];
};

#endif