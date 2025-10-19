//William Whatley
//Objects Header File

//Header guards
#ifndef OBJECTS_H
#define OBJECTS_H

//Other files included
#include "matrix.h"

//libraries included
#include <algorithm>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

class Ray{
    public:
        Ray();
        //assume that the index of reflection is 1.0
        Ray(const Vector3d&, const Vector3d&, int dep = 0, double ir = 1.0);
        Vector3d eye;
        Vector3d dir;
        int ior;
        int depth;
};

class Hit{
    public:
        Hit();
        Hit(double,double,double,double, const Vector3d&, const Vector3d&, const Vector3d&, int);
        void transfer(const double *fil);
        double t{}, alpha{}, beta{}, gamma{};
        Vector3d inter; //determines intersection
        Vector3d norm; //normalized perpendicular to intersection pt
        Vector3d view; //determines the viewing direction
        int raydepth{}; //determines the depth
        double fill[8]{};
};

class Surface{
    public:
        virtual ~Surface();
        //all virtual functions
        virtual bool intersect(const Ray&, double&, double&, Hit&) = 0;
        virtual Surface* clone() = 0;
        virtual void getBounds(Vector3d&, Vector3d&) const = 0;
        double m_fill[8];
        bool isPatch;
};
//Polygon class
class Polygon: public Surface{
    public:
        //constructors/destructors/operators
        Polygon();
        Polygon(int vertices, const double*, bool tri = false, bool pat = false);
        Polygon(const Polygon&);
        virtual ~Polygon();
        Polygon& operator=(const Polygon&);
        Surface* clone();

        //important functions
        bool intersect(const Ray&, double &, double &, Hit &) override;
        void getBounds(Vector3d&, Vector3d&) const override;

        //other functions and member variables
        void pushBackVector(const Vector3d&);
        void addVertex(double x, double y, double z);
        int m_verticies{};
        std::vector<Vector3d> m_vertex;
        bool isTriangle{};
};

//patch class
//child of Polygon
class Patch: public Polygon{
    public:
        //constructor/destructor/assignment operator
        Patch(int, const double*, bool isTri = false, bool pat = true);
        Patch(const Patch& rhs);
        ~Patch();
        Patch& operator=(const Patch&);
        Surface* clone() override;

        //functions
        Vector3d interpolate(const Hit&);
        void addNormal(double, double, double);
        void pushBackNorm(const Vector3d&);
        //member variables
        std::vector<Vector3d> m_normal;
};

//Sphere class
class Sphere: public Surface{
    public:
        //constructor/destructor/operator
        Sphere();
        Sphere(const Sphere&);
        Sphere(double, double, double, double, const double*, bool pat = false);
        ~Sphere();
        Sphere& operator=(const Sphere&);
        Surface* clone();

        //other functions and member variables
        bool intersect(const Ray&, double&, double&, Hit&) override;
        void getBounds(Vector3d&, Vector3d&) const override;
        Vector3d coords;
        double radius;
};

//Node class
class Node{
    public:
    //functions and member variables
        Node(std::vector<Surface*>&, unsigned int, unsigned int);
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


//Viewpoint class
class Viewpoint{
    public:
        //viewpoint constructor
        Viewpoint();

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

//Light class
class Light{
    public:
        Light();
        Light(double, double, double);
        //Vector3d intensity();
        Vector3d coords;

};

#endif