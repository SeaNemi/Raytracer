//William Whatley
//Scene Header file

//header guards
#ifndef SCENE_H
#define SCENE_H

//includes the other file with the objects necessary for the scene
#include "Objects.h"

//libraries included
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

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
        //constructors/operators
        //not  having this assignment operator caused me more pain than you would ever know :(
        Scene();
        Scene(const std::string&);
        Scene(const Scene&);
        
        //destructor and clear, recursiveClear private member function
        ~Scene();
        void clear();
        Scene& operator=(const Scene&);

        //public code
        void setFile(const std::string&);
        void clearFile();
        bool parse();

        //member variables
        Vector3d m_background;
        std::vector<Surface*> m_surfaces;
        std::vector<Light*> m_lights;
        Viewpoint* m_viewpoint;
        Node* m_root;

    private:
        //helper functions
        //for parsing and deletion
        void recursiveClear(Node* node);
        bool parseNFF();
        bool parseObj();

        //parser earclipping and earclipping helpers
        std::vector<Polygon*> earclip(const Polygon*);
        std::vector<Patch*> earclip(const Patch* currPatch);
        bool pointInTriangle(const Vector2d&, const Vector2d&, const Vector2d&, const Vector2d&);
        bool isClockwise(const std::vector<Vector2d>&);
        
        //extra helpers and m_file member variable
        void setBackground(double, double, double);
        void addLight(double, double, double);
        void addSphere(double, double, double, double, const double*);
        std::string m_file;
};
#endif