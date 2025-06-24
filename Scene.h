//William Whatley
//Header file for Scene class

//header guards
#ifndef SCENE_H
#define SCENE_H

//includes the other file with the objects necessary for the scene
#include "Objects.h"

//libraries included
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

//enum used to determine what state the parser is in
enum class State{
    NORMAL,
    VIEWPOINT,
    DETERMINE,
    SPHERE,
    POLYGON,
    PATCH
};

//Node class
//private except only to binary tree
class Node{
    public:
    bool boxIntersect(const Vector3d&, const Vector3d&, const Ray&, double, double);
    void Node::surroundingBox(const Vector3d&, const Vector3d&, const Vector3d&, const Vector3d&, Vector3d&, Vector3d&);
    //member variables
    Node* m_left;
    Node* m_right;
    Surface* m_data;
    Vector3d boxMin;
    Vector3d boxMax;

};

class Scene{
    public:
        Scene();
        Scene(const std::string&);
        Scene(const Scene& rhs);
        ~Scene();
        void setFile(const std::string&);
        void clearFile();
        bool parse();

        //not  having this assignment operator caused me more pain
        //than you would ever know :(
        Scene& operator=(const Scene&);
        void clear();

        //member variables
        Vector3d m_background;
        std::vector<Surface*> m_surfaces;
        std::vector<Light*> m_lights;
        Viewpoint* m_viewpoint;

    private:
        void setBackground(double, double, double);
        void addLight(double, double, double);
        void addSphere(double, double, double, double, const double*);
        bool parseNFF();
        bool parseObj();
        std::vector<Polygon*> earclip(const Polygon*);
        std::vector<Patch*> earclip(const Patch* currPatch);
        bool pointInTriangle(const Vector2d&, const Vector2d&, const Vector2d&, const Vector2d&);
        bool isClockwise(const std::vector<Vector2d>&);

        //member variables
        std::string m_file;

};
#endif