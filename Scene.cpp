//William Whatley
//Scene functions are here

#include "Scene.h"


//Constructor
Scene::Scene(){
    m_file = "NONE";
    m_viewpoint = new Viewpoint();
}
//Overloaded constructor
Scene::Scene(const std::string& file){
    m_file = file;
    m_viewpoint = new Viewpoint();
}

//Destructor
//Cleans all dynamically allocated values
Scene::~Scene(){
    clear();
}

void Scene::clear(){
    //deletes based on whether or not it was created in the first place
    for (Surface* surface : m_surfaces){
        delete surface;
    }


    for (Light *light: m_lights){
        delete light;
    }

    //clears for safety
    m_surfaces.clear();
    m_lights.clear();


    delete m_viewpoint;
    m_viewpoint = nullptr;
}

//sets new file
void Scene::setFile(const std::string& file){
    m_file = file;
}

//clears the current file saved
void Scene::clearFile(){
    m_file = "NONE";
}

//the parser itself
bool Scene::parse(){
    //if no file set up, then return
    if(m_file == "None"){
        return false;
    }

    //opens the io stream to read through the file
    std::ifstream file(m_file);
    std::string line;

    //if the file can't be opened, return
    if(!file.is_open()){
        return false;
    }

    State currState = State::NORMAL;

    //used to determine polygons and fill type
    int verticies;
    int counter = 0;
    double currFill[8];
    Polygon* currPolygon = nullptr;
    Patch* currPatch = nullptr;

    while(std::getline(file, line)){
        std::istringstream curr(line);
        char token;
        std::string word;

        //checks to see what state the currState is in
        switch(currState){
            //if it is normal, that means that no flags are raised
            //thus, parses for beginning letter tokens
            case State::NORMAL:
                curr >> word;

                token = word[0];

                //switch statement determines what the token is
                switch(token){
                    //in case b, three doubles are read in
                    //then added to the background via function call
                    case 'b':
                        double r, g, b;
                        curr >> r >> g >> b;
                        setBackground(r,g,b);
                        break;
                    //if case l, then a light source is being added
                    case 'l':
                        double x, y, z;
                        curr >> x >> y >> z;
                        addLight(x, y, z);
                        break;
                    //if case v, then it is about to be viewpoints
                    //need to account for incoming tokens
                    case 'v':
                        currState = State::VIEWPOINT;
                        break;
                    //if it is f, then a sphere or polygon is coming next
                    //need to account for that
                    //fills m_fill which will then be used for polygons
                    case 'f':
                        currState = State::DETERMINE;
                        for (int i = 0; i < 8; i++){
                            curr >> currFill[i];
                        }
                        break;
                    default:
                        break;
                }
                //end of the NORMAL case
                break;

            //next is the viewpoint case which is concerned with getting the viewpoint
            case State::VIEWPOINT:
                //word is now the new token
                curr >> word;

                //if statement determines which token it is
                //changes m_viewpoint to the new data
                if(word == "from"){
                    double x, y, z;
                    curr >> x >> y >> z;
                    m_viewpoint->setFrom(x, y, z);

                }else if (word == "at"){
                    double x, y, z;
                    curr >> x >> y >> z;
                    m_viewpoint->setAt(x, y, z);

                }else if (word == "up"){
                    double x, y, z;
                    curr >> x >> y >> z;
                    m_viewpoint->setUp(x, y, z);
                }else if (word == "angle"){
                    double g;
                    curr >> g;
                    m_viewpoint->setAngle(g);
                }else if (word == "hither"){
                    double g;
                    curr >> g;
                    m_viewpoint->setHither(g);
                }else if(word == "resolution"){
                    //else if resolution used instead of else to
                    //account for whitespace and comments if influenced
                    int a, b;
                    curr >> a >> b;
                    m_viewpoint->setRes(a, b);

                    //state set back to normal since end of VIEWPOINT stream
                    currState = State::NORMAL;
                }
                break;
            //next is the case of the determine
            //determines if it is a polygon or sphere next
            case (State::DETERMINE):
                curr >> word;
                token = word[0];

                //switches based on the token
                switch(token){
                    case 'f':
                        for (int i = 0; i < 8; i++){
                            curr >> currFill[i];
                        }
                        break;
                    //in case s, reads in the sphere
                    case 's':
                        double x, y, z, r;
                        curr >> x >> y >> z >> r;
                        addSphere(x, y, z, r, currFill);
                        break;
                    
                    //in case p, determines number of veritices
                    //then creates a new polygon and switches state
                    case 'p':
                        counter = 0;
                        if(curr.str()[1] == 'p'){
                            curr >> verticies;
                            //if the second character is a p, means that it is a patch
                            currState = State::PATCH;
                            currPatch = new Patch(verticies, currFill);
                            if(verticies == 3){
                                currPatch->isTriangle = true;
                            }
                        }
                        else{
                            //else, it is just a regular polygon
                            curr >>  verticies;
                            currPolygon = new Polygon(verticies, currFill);
                            currState = State::POLYGON;
                            if(verticies == 3){
                                currPolygon->isTriangle = true;
                            }
                        }
                    default:
                        break;
                }
                break;
            case (State::POLYGON):
                double x, y, z;
                curr >> x >> y >> z;
                currPolygon->addVertex(x, y, z);
                counter++;

                //if all verticies are read in, then the state goes
                //back to normal, and the polygon is then added
                //to the vector of polygons
                if ((counter == verticies)){
                    bool makeTriangles = false;
                    if (verticies == 4) {
                        Vector3d n0 = cross(currPolygon->m_vertex[1] - currPolygon->m_vertex[0], currPolygon->m_vertex[2] - currPolygon->m_vertex[0]);
                        Vector3d n1 = cross(currPolygon->m_vertex[2] - currPolygon->m_vertex[1], currPolygon->m_vertex[3] - currPolygon->m_vertex[1]);
                        Vector3d n2 = cross(currPolygon->m_vertex[3] - currPolygon->m_vertex[2], currPolygon->m_vertex[0] - currPolygon->m_vertex[2]);
                        Vector3d n3 = cross(currPolygon->m_vertex[0] - currPolygon->m_vertex[3], currPolygon->m_vertex[1] - currPolygon->m_vertex[3]);
                        if (((n0.dot(n1) > 0) && ((n0.dot(n2)) > 0) && (n0.dot(n3) > 0)) || (currPolygon->isTriangle)) {
                            makeTriangles = true;
                        }
                    }
                        if (!makeTriangles) {
                            m_surfaces.push_back(new Polygon(*currPolygon));
                        }
                        else{
                            //supports triangle fanning
                            for (int i = 1; i < (verticies - 1); i++){
                                Polygon* ptr = new Polygon(3, currPolygon->m_fill, true);
                                
                                //for normal sides
                                ptr->pushBackVector(currPolygon->m_vertex[0]);
                                ptr->pushBackVector(currPolygon->m_vertex[i]);
                                ptr->pushBackVector(currPolygon->m_vertex[i + 1]);

                                //push back to the m_surfaces
                                m_surfaces.push_back(ptr);
                            }
                        }
                    currState = State::DETERMINE;
                    if(currPolygon){
                        delete currPolygon;
                        currPolygon = nullptr;
                    }
                }

                break;
            
            //case for when a patch is introduced
            case (State::PATCH):
                double a, b, c, d, e, f;
                curr >> a >> b >> c >> d >> e >> f;
                currPatch->addVertex(a, b, c);
                currPatch->addNormal(d, e, f);
                counter++;

                //checks if counter should be reset
                if ((counter == verticies)){
                    bool makeTriangles = false;
                    if (verticies == 4) {
                        Vector3d n0 = cross(currPatch->m_vertex[1] - currPatch->m_vertex[0], currPatch->m_vertex[2] - currPatch->m_vertex[0]);
                        Vector3d n1 = cross(currPatch->m_vertex[2] - currPatch->m_vertex[1], currPatch->m_vertex[3] - currPatch->m_vertex[1]);
                        Vector3d n2 = cross(currPatch->m_vertex[3] - currPatch->m_vertex[2], currPatch->m_vertex[0] - currPatch->m_vertex[2]);
                        Vector3d n3 = cross(currPatch->m_vertex[0] - currPatch->m_vertex[3], currPatch->m_vertex[1] - currPatch->m_vertex[3]);
                        if (((n0.dot(n1) > 0) && ((n0.dot(n2)) > 0) && (n0.dot(n3) > 0)) || (currPatch->isTriangle)) {
                            makeTriangles = true;
                        }
                    }
                    if (!makeTriangles) {
                        m_surfaces.push_back(new Patch(*currPatch));
                    }
                    else{
                        //supports triangle fanning
                        for (int i = 1; i < (verticies - 1); i++){
                            Patch* ptr = new Patch(3, currPatch->m_fill, true);

                            ptr->pushBackVector(currPatch->m_vertex[0]);
                            ptr->pushBackVector(currPatch->m_vertex[i]);
                            ptr->pushBackVector(currPatch->m_vertex[i + 1]);

                            //now does the normalized vectors
                            ptr->pushBackNorm((currPatch->m_normal[0].normalized()));
                            ptr->pushBackNorm(currPatch->m_normal[i].normalized());
                            ptr->pushBackNorm(currPatch->m_normal[i + 1].normalized());

                            //push back to the m_surfaces
                            m_surfaces.push_back(ptr);
                        }
                    }
                    if(currPatch){
                        //delete patch, set to nullptr, continue  
                        delete currPatch;
                        currPatch = nullptr;
                    }
                    currState = State::DETERMINE;
                }
                
                break;
            //default shouldn't be hit but exists for safety
            default:
                currState = State::NORMAL;
                break;
        }
    }

    //the parsing is then FINALLY done
    file.close();
    return true;
}

void Scene::setBackground(double r, double g, double b){
    //vector created then set as the background
    Vector3d vector(r, g, b);
    m_background = vector;
}

//adds a light source to the m_lights vector

void Scene::addLight(double x, double y, double z){
    Light* light = new Light(x, y, z);
    m_lights.push_back(light);
}

//adds a sphere to the m_sphere vector
void Scene::addSphere(double x, double y, double z, double r, const double* fill){
    Sphere* sphere = new Sphere(x, y, z, r, fill);
    m_surfaces.push_back(sphere);
}


//assignment operator
Scene& Scene::operator=(const Scene& rhs){
    if (this != &rhs){
        clear();

        m_file = rhs.m_file;
        m_viewpoint = nullptr;

        for (int i = 0; i < 3; i++){
            m_background[i] = rhs.m_background[i];
        }

        //transfers over the lights
        for(unsigned int i = 0; i < rhs.m_lights.size(); i++){
            Light* light = new Light(rhs.m_lights[i]->coords.x(), rhs.m_lights[i]->coords.y(), rhs.m_lights[i]->coords.z());
            m_lights.push_back(light);
        }


        //transfers over the spheres
        for(const auto& s : rhs.m_surfaces){
            m_surfaces.push_back(s->clone());
        }
        //copies contects of m_viewpoint over
        if(rhs.m_viewpoint){
            m_viewpoint = new Viewpoint();
            m_viewpoint->at = rhs.m_viewpoint->at;
            m_viewpoint->from = rhs.m_viewpoint->from;
            m_viewpoint->up = rhs.m_viewpoint->up;
            m_viewpoint->angle = rhs.m_viewpoint->angle;
            m_viewpoint->hither = rhs.m_viewpoint->hither;

            for (int i = 0; i < 2; i++){
                m_viewpoint->res[i] = rhs.m_viewpoint->res[i];
            }
        }
    }

    return *this;
}