//William Whatley
//Scene functions are here
#include "Scene.h"

//Constructor
Scene::Scene(){
    m_file = "NONE";
    m_root = nullptr;
    m_viewpoint = new Viewpoint();
}
//Overloaded constructor
Scene::Scene(const std::string& file){
    m_file = file;
    m_root = nullptr;
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

    recursiveClear(m_root);
    delete m_viewpoint;
    m_viewpoint = nullptr;
}

//recursiveClear
//uses recursion to completely delete the tree, using postorder traversal
void Scene::recursiveClear(Node* node){
    if (node != nullptr){
        recursiveClear(node->m_left);
        recursiveClear(node->m_right);
        delete node;
    }
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
    if(m_file == "None") return false;
    return (m_file.substr(m_file.find_last_of(".") + 1) == "obj") ? parseObj() : parseNFF();
}

//parseNFF
//parses nff files
bool Scene::parseNFF(){
    //opens the io stream to read through the file
    std::ifstream file(m_file);
    std::string line;

    //if the file can't be opened, return
    if(!file.is_open()) return false;

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
                        for (short i = 0; i < 8; i++){
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
                        for (short i = 0; i < 8; i++){
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
                    //if triangle, skip check
                    if(currPolygon->isTriangle) m_surfaces.push_back(new Polygon(*currPolygon));
                    //else, perform ear clipping algorithm
                    else{
                        std::vector<Polygon*> clipped = earclip(currPolygon);
                        for (unsigned int i = 0; i < clipped.size(); i++) m_surfaces.push_back(clipped[i]);                        
                    }            
                    //now reset the state and delete the ptr
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
                //if all verticies are read in, then the state goes
                //back to normal, and the polygon is then added
                //to the vector of polygons
                if ((counter == verticies)){
                    //if triangle, skip check
                    if(currPatch->isTriangle) m_surfaces.push_back(new Patch(*currPatch));
                    //else, perform ear clipping algorithm
                    else{
                        std::vector<Patch*> clipped = earclip(currPatch);
                        for (unsigned int i = 0; i < clipped.size(); i++) m_surfaces.push_back(clipped[i]);                        
                    }            
                    //now reset the state and delete the ptr
                    currState = State::DETERMINE;
                    if(currPatch){
                        delete currPatch;
                        currPatch = nullptr;
                    }
                }
                break;
            //default shouldn't be hit but exists for safety
            default:
                currState = State::NORMAL;
                break;
        }
    }

    //the parsing is then FINALLY done, then create the root
    file.close();
    m_root = new Node(m_surfaces, 0, m_surfaces.size());
    return true;
}


bool Scene::parseObj(){
    /*
    //opens the io stream to read through the file
    std::ifstream file(m_file);
    std::string line;

    //if the file can't be opened, return
    if(!file.is_open()) return false;

    while(std::getline(file, line)){
        std::istringstream curr(line);
        char token;
        std::string word;







    }
    //the parsing is then FINALLY done
    file.close();
    */
    return true;
}

//pointInTriangle
//used for a point in triangle test
bool Scene::pointInTriangle(const Vector2d& p, const Vector2d& a, const Vector2d& b, const Vector2d& c){
    double alpha = ((b.y() - c.y()) * (p.x() - c.x()) + (c.x() - b.x()) * (p.y() - c.y())) / ((b.y() - c.y()) * (a.x() - c.x()) + (c.x() - b.x()) * (a.y() - c.y()));
    double beta  = ((c.y() - a.y()) * (p.x() - c.x()) + (a.x() - c.x()) * (p.y() - c.y())) / ((b.y() - c.y()) * (a.x() - c.x()) + (c.x() - b.x()) * (a.y() - c.y()));
    double gamma = 1.0 - alpha - beta;
    return alpha >= 0 && beta >= 0 && gamma >= 0;
}

//isClockwise
//determines orientation to ensure clockwise position
bool Scene::isClockwise(const std::vector<Vector2d>& verts){
    double sum = 0.0;
    for (unsigned int i = 0; i < verts.size(); ++i) {
        const Vector2d& a = verts[i];
        const Vector2d& b = verts[(i + 1) % verts.size()];
        sum += (b.x() - a.x()) * (b.y() + a.y());
    }
    return sum > 0.0; // CW if positive
}


//earclip
//performs triangluation using the earclip algorithm
std::vector<Polygon*> Scene::earclip(const Polygon* currPoly){
    //result initalized, vertexes copied over from the Polygon sent in
    //size check performed to prevent illegal shaping
    std::vector<Polygon*> result;
    std::vector<Vector3d> vertexes = currPoly->m_vertex;
    if (vertexes.size() < 3) return result;

    //axis project determined based on the polygon's normal
    int projectDir = 2;
    Vector3d normal(0, 0, 0);
    for (unsigned int i = 0; i < vertexes.size(); ++i) {
        const Vector3d& curr = vertexes[i];
        const Vector3d& next = vertexes[(i + 1) % vertexes.size()];
        normal += curr.cross(next);
    }
    normal.normalize();

    //choose correct projection axis
    //discard axis with largest magnitude in the normal vector
    if (fabs(normal[0]) > fabs(normal[1]) && fabs(normal[0]) > fabs(normal[2])) projectDir = 0;
    else if (fabs(normal[1]) > fabs(normal[2])) projectDir = 1;

    //3D vertices projected onto a 2D plane
    std::vector<Vector2d> vexes2d;
    for (const auto& v : vertexes)
        vexes2d.emplace_back(v, projectDir);

    //ensure vertices are in counter-clockwise order since required by ear clipping
    if (isClockwise(vexes2d)) {
        std::reverse(vertexes.begin(), vertexes.end());
        std::reverse(vexes2d.begin(), vexes2d.end());
    }

    //ears are clipped until fewer than 3 vertices remain
    while (vertexes.size() >= 3) {
        bool earFound = false;
        int n = vertexes.size();

        //for loop goes through remaining verticies
        for (int i = 0; i < n; ++i) {
            double pIndex = (i + n - 1) % n;
            double nIndex = (i + 1) % n;

            //get prev, curr, next verticies
            const Vector2d& prev = vexes2d[pIndex];
            const Vector2d& curr = vexes2d[i];
            const Vector2d& next = vexes2d[nIndex];

            //check to see if the corner is convex
            double crossZ = (curr - prev).cross(next - curr);
            if (crossZ <= 0) continue;  // Not convex, skip

            //check to see if any other vertex lies inside the triangle
            //done to prevent inappropriate clipping
            bool pointInside = false;
            for (int j = 0; j < n; ++j) {
                if (j == pIndex || j == i || j == nIndex) continue;
                if (pointInTriangle(vexes2d[j], prev, curr, next)) {
                    pointInside = true;
                    break;
                }
            }
            if (pointInside) continue;

            //if no other vertex is inside and the corner is convex, we found an ear
            Polygon* tri = new Polygon(3, currPoly->m_fill, true);
            tri->pushBackVector(vertexes[pIndex]);
            tri->pushBackVector(vertexes[i]);
            tri->pushBackVector(vertexes[nIndex]);
            result.push_back(tri);

            //remove the ear tip from the polygon
            vertexes.erase(vertexes.begin() + i);
            vexes2d.erase(vexes2d.begin() + i);
            earFound = true;
            //break and restart search with new verticies removed
            break;
        }

        //if we have no ear, break to prevent infinite loop
        if (!earFound) break;
    }

    //forced fan triangulation used as a last resort
    //forces triangulations 
    if (vertexes.size() >= 3) {
        //used to determine orientation
        Vector3d sumNormal(0, 0, 0);
        for (unsigned int i = 0; i < vertexes.size(); ++i) {
            const Vector3d& v0 = vertexes[i];
            const Vector3d& v1 = vertexes[(i + 1) % vertexes.size()];
            sumNormal += v0.cross(v1);
        }
        sumNormal.normalize();

        //flip vertex order if the polygon is wound downward
        if (sumNormal.z() < 0) std::reverse(vertexes.begin(), vertexes.end());

        //create forced triangle fan from the first vertex
        for (unsigned int i = 1; i + 1 < vertexes.size(); ++i) {
            Polygon* tri = new Polygon(3, currPoly->m_fill, true);
            tri->pushBackVector(vertexes[0]);
            tri->pushBackVector(vertexes[i]);
            tri->pushBackVector(vertexes[i + 1]);
            result.push_back(tri);
        }
    }
    return result;
}

//earclip overload for patch
//performs triangluation using the earclip algorithm for smooth shaded patches
std::vector<Patch*> Scene::earclip(const Patch* currPatch) {
    //result initalized, vertexes and normals copied over from the Patch sent in
    //size check performed to prevent illegal shaping
    std::vector<Patch*> result;
    std::vector<Vector3d> vertexes = currPatch->m_vertex;
    std::vector<Vector3d> normals = currPatch->m_normal;
    if (vertexes.size() < 3) return result;

    //axis project determined based on the polygon's normal
    int projectDir = 2;
    Vector3d normal(0, 0, 0);
    for (unsigned int i = 0; i < vertexes.size(); ++i) {
        const Vector3d& curr = vertexes[i];
        const Vector3d& next = vertexes[(i + 1) % vertexes.size()];
        normal += curr.cross(next);
    }
    normal.normalize();

    //choose correct projection axis
    //discard axis with largest magnitude in the normal vector
    if (fabs(normal[0]) > fabs(normal[1]) && fabs(normal[0]) > fabs(normal[2])) projectDir = 0;
    else if (fabs(normal[1]) > fabs(normal[2])) projectDir = 1;

    //3D vertices projected onto a 2D plane
    std::vector<Vector2d> vexes2d;
    for (const auto& v : vertexes)
        vexes2d.emplace_back(v, projectDir);

    //ensure vertices are in counter-clockwise order since required by ear clipping
    if (isClockwise(vexes2d)) {
        std::reverse(vertexes.begin(), vertexes.end());
        std::reverse(normals.begin(), normals.end());
        std::reverse(vexes2d.begin(), vexes2d.end());
    }

    //ears are clipped until fewer than 3 vertices remain
    while (vertexes.size() >= 3) {
        bool earFound = false;
        int n = vertexes.size();

        //for loop goes through remaining verticies
        for (int i = 0; i < n; ++i) {
            int pIndex = (i + n - 1) % n;
            int nIndex = (i + 1) % n;

            //get prev, curr, next verticies
            const Vector2d& prev = vexes2d[pIndex];
            const Vector2d& curr = vexes2d[i];
            const Vector2d& next = vexes2d[nIndex];

            //check to see if the corner is convex
            double crossZ = (curr - prev).cross(next - curr);
            if (crossZ <= 0) continue;  // Not convex, skip

            //check to see if any other vertex lies inside the triangle
            //done to prevent inappropriate clipping
            bool pointInside = false;
            for (int j = 0; j < n; ++j) {
                if (j == pIndex || j == i || j == nIndex) continue;
                if (pointInTriangle(vexes2d[j], prev, curr, next)) {
                    pointInside = true;
                    break;
                }
            }
            if (pointInside) continue;

            //if no other vertex is inside and the corner is convex, we found an ear
            Patch* tri = new Patch(3, currPatch->m_fill, true);
            tri->pushBackVector(vertexes[pIndex]);
            tri->pushBackVector(vertexes[i]);
            tri->pushBackVector(vertexes[nIndex]);

            tri->pushBackNorm(normals[pIndex]);
            tri->pushBackNorm(normals[i]);
            tri->pushBackNorm(normals[nIndex]);

            result.push_back(tri);

            //remove the ear tip from the polygon
            vertexes.erase(vertexes.begin() + i);
            normals.erase(normals.begin() + i);
            vexes2d.erase(vexes2d.begin() + i);
            earFound = true;
            //break and restart search with new verticies removed
            break;
        }
        //if we have no ear, break to prevent infinite loop
        if (!earFound) break;
    }

    //forced fan triangulation used as a last resort
    //forces triangulations
    if (vertexes.size() >= 3) {
        //used to determine orientation
        Vector3d sumNormal(0, 0, 0);
        for (unsigned int i = 0; i < vertexes.size(); ++i) {
            const Vector3d& v0 = vertexes[i];
            const Vector3d& v1 = vertexes[(i + 1) % vertexes.size()];
            sumNormal += v0.cross(v1);
        }
        sumNormal.normalize();

        //flip vertex order if the polygon is wound downward
        if (sumNormal.z() < 0) {
            std::reverse(vertexes.begin(), vertexes.end());
            std::reverse(normals.begin(), normals.end());
        }

        //create forced triangle fan from the first vertex
        for (unsigned int i = 1; i + 1 < vertexes.size(); ++i) {
            Patch* tri = new Patch(3, currPatch->m_fill, true);
            tri->pushBackVector(vertexes[0]);
            tri->pushBackVector(vertexes[i]);
            tri->pushBackVector(vertexes[i + 1]);
            tri->pushBackNorm(normals[0]);
            tri->pushBackNorm(normals[i]);
            tri->pushBackNorm(normals[i + 1]);
            result.push_back(tri);
        }
    }
    return result;
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
        m_root = nullptr;

        for (short i = 0; i < 3; i++){
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

            for (short i = 0; i < 2; i++){
                m_viewpoint->res[i] = rhs.m_viewpoint->res[i];
            }
        }
        m_root = new Node(m_surfaces, 0, m_surfaces.size());
    }

    return *this;
}