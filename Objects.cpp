//William Whatley
//Objects source file

#include "Objects.h"

//ray class
Ray::Ray(Vector3d e, Vector3d d, int dep){
    eye = e;
    dir = d;
    depth = dep;
}

//added this so that way the Raytracer can actually read it
Hit::Hit(double to,double a,double b,double g, Vector3d inte, Vector3d nor, Vector3d vie,  int depth){
    t = to;
    alpha = a;
    beta = b;
    gamma = g;
    inter = inte;
    view = vie;
    norm = nor;
    raydepth = depth;
}

void Hit::transfer(double fil[8]){
    for(int i = 0; i < 8; i++){
        fill[i] = fil[i];
    }
}

//Surface destructor
//needed for destruction of child classes
Surface::~Surface(){

}

//Constructor for Polygon
Polygon::Polygon(int verticies, double fill[8], bool tri, bool pat){
    m_verticies = verticies;
    for (int i = 0; i < 8; i++){
        m_fill[i] = fill[i];
    }

    isTriangle = tri;

    isPatch = pat;
}

//deletes the vertex sides
Polygon::~Polygon(){
    clear();
}

void Polygon::clear(){
    for (unsigned int i = 0; i < m_vertex.size(); i++){
        delete m_vertex[i];
        m_vertex[i] = nullptr;
    }
}
//adds a vertex side to the polygon
void Polygon::addVertex(double x, double y, double z){
    Vector3d* vertex = new Vector3d(x, y, z);
    m_vertex.push_back(vertex);
}

Polygon::Polygon(const Polygon& rhs){
    m_verticies = rhs.m_verticies;
    //sets up everything
    for (int i = 0; i < 8; i++){
        m_fill[i] = rhs.m_fill[i];
    }
    
    //for loop pushes back the vectors
    for (const Vector3d* vertex : rhs.m_vertex){
        if(vertex == nullptr){
            m_vertex.push_back(nullptr);
        }else{
            addVertex(vertex->x(), vertex->y(), vertex->z());
        }
    }
    
    isPatch = rhs.isPatch;
    isTriangle = rhs.isTriangle;
}

//assignment operator
Polygon& Polygon::operator=(const Polygon& rhs){
    //checks for self assignment
    if(this != &rhs){
        clear();
        m_verticies = rhs.m_verticies;
        //sets up everything
        for (int i = 0; i < 8; i++){
            m_fill[i] = rhs.m_fill[i];
        }
        
        //for loop pushes back the vectors
        for (const Vector3d* vertex : rhs.m_vertex){
            if(vertex == nullptr){
                m_vertex.push_back(nullptr);
            }else{
                addVertex(vertex->x(), vertex->y(), vertex->z());
            }
        }
    }
    return *this;
}

//Polygon intersect
//determines if intersects with polygon
bool Polygon::intersect(Ray ray, double &min, double &max, Hit &hr){
    if(isTriangle){
        //first matrix set up is matrix a
        Matrix3d matrixA;

        //the columns are then set based on the book
        matrixA.col(0, (*m_vertex[0] - *m_vertex[1]));
        matrixA.col(1, (*m_vertex[0] - *m_vertex[2]));
        matrixA.col(2, ray.dir);

        //determinant taken
        double detA = matrixA.determinant();

        if (detA == 0){
            return false;
        }
        
        //first value to be solved for is t
        //need to create matrix and then divide
        Matrix3d matrixT;
        matrixT.col(0, (*m_vertex[0] - *m_vertex[1]));
        matrixT.col(1, (*m_vertex[0] - *m_vertex[2]));
        matrixT.col(2, (*m_vertex[0] - ray.eye));
        double t = matrixT.determinant() / detA;

        //that is then compared to the limits
        if ((t < min) || (t > max)){
            return false;
        }

        //next need to solve for upsilon
        //next to be solved for is matrixUpsilon
        Matrix3d matrixUpsilon;
        matrixUpsilon.col(0, (*m_vertex[0] - *m_vertex[1]));
        matrixUpsilon.col(1, (*m_vertex[0] - ray.eye));
        matrixUpsilon.col(2, ray.dir);

        double upsilon = matrixUpsilon.determinant() / detA;

        //then we check if it is valid
        if ((upsilon < 0) || (upsilon > 1)){
            return false;
        }

        //finally we need to solve for beta
        Matrix3d matrixBeta;
        matrixBeta.col(0, (*m_vertex[0] - ray.eye));
        matrixBeta.col(1, (*m_vertex[0] - *m_vertex[2]));
        matrixBeta.col(2, ray.dir);

        double beta = matrixBeta.determinant() / detA;

        //we then compare beta and see if it is valid
        //if it is, that means we have indeed hit a triangle
        if ((beta < 0) || (beta + upsilon > 1)){
            return false;
        }

        //taken from professsor's code
        max = t;
        hr.t = t; 
        hr.inter = ray.eye + (t * ray.dir);
        hr.norm = (*m_vertex[0] - *m_vertex[1]).cross((*m_vertex[0] - *m_vertex[2]));
        hr.norm.normalize();
        hr.alpha = 1.0 - beta - upsilon;
        hr.beta = beta;
        hr.gamma = upsilon;
        
        return true;
    }else{
        int projectDir;
        Vector3d n = (*m_vertex[1]- *m_vertex[0]).cross(*m_vertex[2] - *m_vertex[0]);
        n.normalize();

        double t = -(ray.eye.dot(n) - m_vertex[0]->dot(n)) / (ray.dir.dot(n));
        if (t < min || t > max) return false;

        Vector3d p;
        p = ray.eye + t * ray.dir;

        if (fabs(n[0]) > fabs(n[1]) && fabs(n[0]) > fabs(n[2])) {
            projectDir = 0;
        } else if (fabs(n[1]) > fabs(n[2])) {
            projectDir = 1;
        } else {
            projectDir = 2;
        }
        
        //creates two projected vector
        Vector2d p2(p, projectDir);
        Vector2d bbMin(*m_vertex[0], projectDir);

        Vector2d bbMax = bbMin;
        for (unsigned int i=1; i<m_vertex.size(); i++) {
            Vector2d v(*m_vertex[i], projectDir);
            if (v[0] < bbMin[0]) bbMin[0] = v[0];
            if (v[0] > bbMax[0]) bbMax[0] = v[0];
            if (v[1] < bbMin[1]) bbMin[1] = v[1];
            if (v[1] > bbMax[1]) bbMax[1] = v[1];
        }
        
        if (p2[0] < bbMin[0]) return false;
        if (p2[1] < bbMin[1]) return false;
        if (p2[0] > bbMax[0]) return false;
        if (p2[1] > bbMax[1]) return false;

        Vector2d dir(sqrt(2), sqrt(2));
        int count = 0;
        for (unsigned int i=0; i<m_vertex.size(); i++) {
            //two projected vectors are created
            Vector2d a(*m_vertex[i], projectDir);
            Vector2d b(*m_vertex[(i+1) % m_vertex.size()], projectDir);
            Vector2d ab = b-a;
            Vector2d ap2(a-p2);
            double t2 = ap2.cross(ab / dir.cross(ab));
            if (t2 < 0.0) continue;
            double alpha = ap2.cross(dir / dir.cross(ab));
            if (alpha > 0.0 && alpha < 1.0) count++;
        }

        if (count % 2 == 0) {
            return false;
        }

        hr.t = t;
        hr.inter = p;
        hr.norm = n;
        return true;
    }
}

int Polygon::getSize(){
    return m_verticies;
}

//pushes back a vertex already in eigen form
void Polygon::pushBackVector(Vector3d* push){
    m_vertex.push_back(push);
}

Surface* Polygon::clone(){
    return new Polygon(*this);
}

//creates a sphere
Sphere::Sphere(double x, double y, double z, double r, double fill[8], bool pat){
    Vector3d vector(x, y, z);
    coords = vector;
    radius = r;
    for (int i = 0; i < 8; i++){
        m_fill[i] = fill[i];
    }

    isPatch = pat;
}

Sphere::Sphere(const Sphere& rhs){
    Vector3d vector(rhs.coords);
    coords = vector;
    radius = rhs.radius;
    for (int i = 0; i < 8; i++){
        m_fill[i] = rhs.m_fill[i];
    }
    isPatch = rhs.isPatch;
}

Surface* Sphere::clone(){
    return new Sphere(*this);
}

//raySph
//checks to see if a sphere is hit
bool Sphere::intersect(Ray ray, double &min, double &max, Hit &hr){
    //solving for determinant
    //a is direction dot direction

    Vector3d vector = ray.dir - coords;
    Vector3d direction = ray.dir;

    double a = direction.squaredNorm();
    //b is direction dot eye - coords
    double b = direction.dot(vector);

    //c is eye - coords squared, minus radius squared
    double c = vector.squaredNorm()- (radius * radius);

    //first need to solve for discriminant to see if valid
    double discriminant = (b * b) - (a * c);

    //if discriminant is zero, then return false, else a sphere is hit
    if (discriminant < 0){
        return false;
    }
    else{
        //else solve using quadratic and return true
        //changes the min and max
        double root1, root2;

        //roots calculated
        double ans = sqrt(discriminant);
        root1 = ((-1 * b) + ans)/(a);
        root2 = ((-1 * b) - ans)/(a);

        //choose appopriate root
        double t = ((root1 < 0) || (root2 > 0 && (root2 < root1))) ? root1 : root2;

        //check that t is in the valid range
        if((t < min) || (t > max)){
            return false;
        }

        max = t;
        hr.t = t;
        return true;
    }
}


//creates a light
Light::Light(double x, double y, double z){
    Vector3d vector(x, y, z);
    coords = vector; 
}

//All of the setters for viewpoint

void Viewpoint::setAt(double x, double y, double z){
    Vector3d vector(x, y, z);
    at = vector; 
}

void Viewpoint::setFrom(double x, double y, double z){
    Vector3d vector(x, y, z);
    from = vector; 
}
void Viewpoint::setUp(double x, double y, double z){
    Vector3d vector(x, y, z);
    up = vector;
}
void Viewpoint::setAngle(double a){
    angle = a;
}
void Viewpoint::setHither(double g){
    hither = g;
}
void Viewpoint::setRes(int x, int y){
    res[0] = x;
    res[1] = y;
}

//constructor
Patch::Patch(int verticies, double fill[8], bool isTri, bool pat){
    m_verticies = verticies;
    for (int i = 0; i < 8; i++){
        m_fill[i] = fill[i];
    }

    isTriangle = isTri;

    isPatch = pat;
}

//copy constructor, unsure if necessary but safe to have
Patch::Patch(const Patch& rhs) : Polygon(rhs){
    for (const auto& norm: rhs.m_normal){
        m_normal.push_back(new Vector3d(*norm));
    }
}

//destructor
Patch::~Patch(){
    for (const auto& norm: m_normal){
        delete norm;
    }
}

//adds a normal vertex to the Patch
void Patch::addNormal(double x, double y, double z){
    Vector3d* vertex = new Vector3d(x, y, z);
    m_normal.push_back(vertex);
}

//push back norma eigenvector
void Patch::pushBackNorm(Vector3d* a){
    m_normal.push_back(a);
}

Surface* Patch::clone(){
    return new Patch(*this);
}

//interpolate
//smoothly shades the patches
Vector3d Patch::interpolate(Hit hit){
    //Vector3d normal = (1 - hit.beta - hit.gamma) * *m_normal[0] + hit.beta * *m_normal[1] + hit.gamma * *m_normal[2];
    Vector3d vec;
    vec = (1 - hit.beta - hit.gamma) * *m_normal[0] + hit.beta * *m_normal[1] + hit.gamma * *m_normal[2];
    return vec;
}