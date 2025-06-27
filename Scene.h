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



class Scene{
    public:
        Scene();
        Scene(const std::string&);
        void treeCount(Node* node, int& data);
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
        Node* m_root;

    private:
        void setBackground(double, double, double);
        void addLight(double, double, double);
        void addSphere(double, double, double, double, const double*);
        bool parseNFF();
        bool parseObj();
        void recursiveClear(Node* node);
        std::vector<Polygon*> earclip(const Polygon*);
        std::vector<Patch*> earclip(const Patch* currPatch);
        bool pointInTriangle(const Vector2d&, const Vector2d&, const Vector2d&, const Vector2d&);
        bool isClockwise(const std::vector<Vector2d>&);

        //member variables
        std::string m_file;

};
#endif