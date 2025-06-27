//William Whatley
//Objects Header File

//Header guards
#ifndef OBJECTS_H
#define OBJECTS_H

//What's included
#include <iostream>
#include "matrix.h"
#include <algorithm>
#include <string>
#include <random>
#include <limits>
#include <vector>

class Ray{
    public:
        Ray();
        Ray(const Vector3d&, const Vector3d&, int dep = 0);
        Vector3d eye;
        Vector3d dir;
        int depth;
};

class Hit{
    public:
        Hit(){}
        Hit(double,double,double,double, const Vector3d&, const Vector3d&, const Vector3d&, int);
        double t, alpha, beta, gamma;
        Vector3d inter; //determines intersection
        Vector3d norm; //normalized perpendicular to intersection pt
        Vector3d view; //determines the viewing direction
        int raydepth; //determines the depth
        void transfer(const double *fil);
        double fill[8];
};

class Surface{
    public:
        virtual ~Surface();
        virtual bool intersect(const Ray&, double&, double&, Hit&) = 0;
        virtual Surface* clone() = 0;
        virtual void getBounds(Vector3d&, Vector3d&) const = 0;
        double m_fill[8];
        bool isPatch;
};
//Polygon class
class Polygon: public Surface{
    public:
        Polygon(){}
        Polygon(int vertices, const double*, bool tri = false, bool pat = false);
        Polygon(const Polygon&);
        virtual ~Polygon();
        int getSize();
        void addVertex(double x, double y, double z);
        void pushBackVector(const Vector3d&);
        Polygon& operator=(const Polygon&);
        bool intersect(const Ray&, double &, double &, Hit &) override;
        void getBounds(Vector3d&, Vector3d&) const override;
        int m_verticies;
        std::vector<Vector3d> m_vertex;

        Surface* clone();
        bool isTriangle;
};

//patch class
//child of Polygon
class Patch: public Polygon{
    public:
        //functions
        Patch(int, const double*, bool isTri = false, bool pat = true);
        Patch(const Patch& rhs);
        ~Patch();
        void addNormal(double, double, double);
        void pushBackNorm(const Vector3d&);
        //member variables
        std::vector<Vector3d> m_normal;
        Vector3d interpolate(const Hit&);
        Surface* clone() override;
};

//Sphere class
class Sphere: public Surface{
    public:
        //public
        Sphere(){}
        Sphere(const Sphere&);
        Sphere(double, double, double, double, const double*, bool pat = false);
        Sphere& operator=(const Sphere&);
        bool intersect(const Ray&, double&, double&, Hit&) override;
        void getBounds(Vector3d&, Vector3d&) const override;
        Surface* clone();

        //member variables
        Vector3d coords;
        double radius;
};


//Light class
class Light{
    public:
        Light();
        Light(double, double, double);
        //Vector3d intensity();
        Vector3d coords;

};

//Viewpoint class
class Viewpoint{
    public:
        //viewpoint constructor has nothing
        Viewpoint(){}
        Viewpoint(const Viewpoint&);
        Viewpoint& operator=(const Viewpoint&);

        //all setters
        void setAt(double, double, double);
        void setFrom(double, double, double);
        void setUp(double, double, double);
        void setAngle(double);
        void setHither(double);
        void setRes(int, int);


        //doubles for the viewpoint
        Vector3d at;
        Vector3d from;
        Vector3d up;
        double angle;
        double hither;
        double res[2];
};

//Node class
//private except only to binary tree
class Node{
    public:
    Node(std::vector<Surface*>&, unsigned int, unsigned int);
    ~Node();
    bool boxIntersect(const Vector3d&, const Vector3d&, const Ray&, double, double);
    void surroundingBox(const Vector3d&, const Vector3d&, const Vector3d&, const Vector3d&, Vector3d&, Vector3d&);
    bool compare(Surface*, Surface*, short);
    unsigned int nelement(std::vector<Surface*>&, unsigned int, unsigned int, short);
    void quickselect(std::vector<Surface*>&, unsigned int, unsigned int, unsigned int, short);
    //member variables
    Node* m_left;
    Node* m_right;
    Surface* m_data;
    Vector3d m_min;
    Vector3d m_max;

};

#endif