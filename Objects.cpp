//William Whatley
//Objects source file
#include "Objects.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// RAY CLASS /////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//ray constructors
Ray::Ray(){}
Ray::Ray(const Vector3d& e, const Vector3d& d, int dep){
    eye = e;
    dir = d;
    depth = dep;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// HIT CLASS /////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//added this so that way the Raytracer can actually read it
Hit::Hit(){}
Hit::Hit(double to,double a,double b,double g, const Vector3d& inte, const Vector3d& nor, const Vector3d& vie,  int depth){
    t = to;
    alpha = a;
    beta = b;
    gamma = g;
    inter = inte;
    view = vie;
    norm = nor;
    raydepth = depth;
}
//transfer used to move stuff over
void Hit::transfer(const double *fil){
    for(short i = 0; i < 8; i++) fill[i] = fil[i];
    
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// SURFACE CLASS /////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//Surface destructor
//needed for destruction of child classes
Surface::~Surface(){}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// POLYGON CLASS /////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//Constructors for Polygon
Polygon::Polygon(){}
Polygon::Polygon(int verticies, const double *fill, bool tri, bool pat){
    m_verticies = verticies;
    for (short i = 0; i < 8; i++) m_fill[i] = fill[i];
    isTriangle = tri;
    isPatch = pat;
}
Polygon::Polygon(const Polygon& rhs){
    m_verticies = rhs.m_verticies;
    //sets up everything
    for (short i = 0; i < 8; i++) m_fill[i] = rhs.m_fill[i];
    //for loop pushes back the vectors
    for (const Vector3d& vertex : rhs.m_vertex) addVertex(vertex.x(), vertex.y(), vertex.z());
    isPatch = rhs.isPatch;
    isTriangle = rhs.isTriangle;
}

//destructor deletes the vertex sides
Polygon::~Polygon(){}

//assignment operator
Polygon& Polygon::operator=(const Polygon& rhs){
    //checks for self assignment
    if(this != &rhs){
        m_verticies = rhs.m_verticies;
        //sets up everything
        for (short i = 0; i < 8; i++) m_fill[i] = rhs.m_fill[i];
        //for loop pushes back the vectors
        for (const Vector3d& vertex : rhs.m_vertex) addVertex(vertex.x(), vertex.y(), vertex.z());
    }
    return *this;
}
//clone
Surface* Polygon::clone(){
    return new Polygon(*this);
}

//Polygon intersect
//determines if intersects with polygon
bool Polygon::intersect(const Ray& ray, double &min, double &max, Hit &hr){
    //if statement determines if triangle or not
    //both use different algorithms
    if(isTriangle){
        //first matrix set up is matrix a
        Matrix3d matrixA;
        //the columns are then set based on the book
        matrixA.col(0, (m_vertex[0] - m_vertex[1]));
        matrixA.col(1, (m_vertex[0] - m_vertex[2]));
        matrixA.col(2, ray.dir);

        //determinant taken
        double detA = matrixA.determinant();
        if (detA == 0)return false;
        
        //first value to be solved for is t
        //need to create matrix and then divide
        Matrix3d matrixT;
        matrixT.col(0, (m_vertex[0] - m_vertex[1]));
        matrixT.col(1, (m_vertex[0] - m_vertex[2]));
        matrixT.col(2, (m_vertex[0] - ray.eye));
        double t = matrixT.determinant() / detA;

        //that is then compared to the limits
        if ((t < min) || (t > max)) return false;

        //next need to solve for upsilon
        //next to be solved for is matrixUpsilon
        Matrix3d matrixUpsilon;
        matrixUpsilon.col(0, (m_vertex[0] - m_vertex[1]));
        matrixUpsilon.col(1, (m_vertex[0] - ray.eye));
        matrixUpsilon.col(2, ray.dir);
        double upsilon = matrixUpsilon.determinant() / detA;

        //then we check if it is valid
        if ((upsilon < 0) || (upsilon > 1))return false;

        //finally we need to solve for beta
        Matrix3d matrixBeta;
        matrixBeta.col(0, (m_vertex[0] - ray.eye));
        matrixBeta.col(1, (m_vertex[0] - m_vertex[2]));
        matrixBeta.col(2, ray.dir);
        double beta = matrixBeta.determinant() / detA;

        //we then compare beta and see if it is valid
        //if it is, that means we have indeed hit a triangle
        if ((beta < 0) || (beta + upsilon > 1))return false;

        //update the information since true
        max = t;
        hr.t = t; 
        hr.inter = ray.eye + (t * ray.dir);
        hr.norm = (m_vertex[0] - m_vertex[1]).cross((m_vertex[0] - m_vertex[2]));
        hr.norm.normalize();
        hr.alpha = 1.0 - beta - upsilon;
        hr.beta = beta;
        hr.gamma = upsilon;
        return true;
    }else{
        //initializes projectDirection and creates n
        int projectDir;
        Vector3d n = (m_vertex[1]- m_vertex[0]).cross(m_vertex[2] - m_vertex[0]);
        n.normalize();
        double t = -(ray.eye.dot(n) - m_vertex[0].dot(n)) / (ray.dir.dot(n));
        if (t < min || t > max) return false;
        Vector3d p= ray.eye + t*ray.dir;
        
        //determines which axis to project onto
        //choose largest axis for stability
        if (std::fabs(n[0]) > std::fabs(n[1]) && std::fabs(n[0]) > std::fabs(n[2])) projectDir = 0;
        else if (std::fabs(n[1]) > std::fabs(n[2])) projectDir = 1;
        else projectDir = 2;
    
        //creates two projected vector
        Vector2d p2(p, projectDir);
        Vector2d bbMin(m_vertex[0], projectDir);

        Vector2d bbMax = bbMin;
        //computes the bounding box for the 2d
        for (unsigned int i=1; i<m_vertex.size(); i++) {
            Vector2d v(m_vertex[i], projectDir);
            if (v[0] < bbMin[0]) bbMin[0] = v[0];
            if (v[0] > bbMax[0]) bbMax[0] = v[0];
            if (v[1] < bbMin[1]) bbMin[1] = v[1];
            if (v[1] > bbMax[1]) bbMax[1] = v[1];
        }
        //reject based on bounding box
        if (p2[0] < bbMin[0]) return false;
        if (p2[1] < bbMin[1]) return false;
        if (p2[0] > bbMax[0]) return false;
        if (p2[1] > bbMax[1]) return false;
        Vector2d dir(std::sqrt(2), std::sqrt(2));
        int count = 0;

        //for loop goes through each edge of the polygon
        for (unsigned int i=0; i<m_vertex.size(); i++) {
            Vector2d a(m_vertex[i], projectDir);
            Vector2d b(m_vertex[(i+1) % m_vertex.size()], projectDir);
            Vector2d ab = b-a;
            double denom = dir.cross(ab);
            Vector2d ap2(a-p2);
            double t2 = ap2.cross(ab/denom);
            if (t2 < 0.0) continue;
            double alpha = ap2.cross(dir / denom);
            if (alpha > 0.0 && alpha < 1.0) count++;
        }

        //if count is even, means that polygon exited and entered equal amount of times, thus point outside
        if (count % 2 == 0) return false;
        //update info if true
        hr.t = t;
        hr.inter = p;
        hr.norm = n;
        return true;
    }
}

//getBounds
//based on veritices of the polygon
void Polygon::getBounds(Vector3d& min, Vector3d& max) const{
    //start with the min and max vertitices being the first in the vector
    //then go through and figure out the min and max
    min = max = m_vertex[0];
    //using ptr logic for more effecitveness
    for (const Vector3d* ptr = &m_vertex[0]; ptr != &m_vertex[0] + m_vertex.size(); ++ptr) {
        for (short j = 0; j < 3; ++j) {
            if ((*ptr)[j] < min[j]) min[j] = (*ptr)[j];
            if ((*ptr)[j] > max[j]) max[j] = (*ptr)[j];
        }
    }
}

//pushes back a vertex already in eigen form
void Polygon::pushBackVector(const Vector3d& push){
    m_vertex.push_back(push);
}

//adds a vertex side to the polygon
void Polygon::addVertex(double x, double y, double z){
    Vector3d vertex(x, y, z);
    m_vertex.push_back(vertex);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// PATCH CLASS ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//constructor/destructor/assignment operator
Patch::Patch(int verticies, const double* fill, bool isTri, bool pat){
    m_verticies = verticies;
    for (int i = 0; i < 8; i++) m_fill[i] = fill[i];
    isTriangle = isTri;
    isPatch = pat;
}
Patch::Patch(const Patch& rhs) : Polygon(rhs){
    for (const auto& norm: rhs.m_normal) m_normal.push_back(Vector3d(norm));
}
Patch::~Patch(){}
Patch& Patch::operator=(const Patch& rhs){
    if(this != &rhs){
        Polygon::operator=(rhs);
        for(unsigned int i = 0; i < rhs.m_normal.size(); i++) m_normal.push_back(rhs.m_normal[i]);
    }
    return *this;
}
Surface* Patch::clone(){
    return new Patch(*this);
}

//interpolate
//smoothly shades the patches
Vector3d Patch::interpolate(const Hit& hit){
    //Vector3d normal = (1 - hit.beta - hit.gamma) * *m_normal[0] + hit.beta * *m_normal[1] + hit.gamma * *m_normal[2];
    Vector3d vec;
    vec = (1 - hit.beta - hit.gamma) * m_normal[0] + hit.beta * m_normal[1] + hit.gamma * m_normal[2];
    return vec;
}

//adds a normal vertex to the Patch
void Patch::addNormal(double x, double y, double z){
    Vector3d vertex(x, y, z);
    m_normal.push_back(vertex);
}

//push back norma eigenvector
void Patch::pushBackNorm(const Vector3d& a){
    m_normal.push_back(a);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// SPHERE CLASS //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//constructors/destructors/operators
Sphere::Sphere(){}
Sphere::Sphere(double x, double y, double z, double r, const double* fill, bool pat){
    Vector3d vector(x, y, z);
    coords = vector;
    radius = r;
    for (short i = 0; i < 8; i++)m_fill[i] = fill[i];
    isPatch = pat;
}
Sphere::Sphere(const Sphere& rhs){
    Vector3d vector(rhs.coords);
    coords = vector;
    radius = rhs.radius;
    for (short i = 0; i < 8; i++) m_fill[i] = rhs.m_fill[i];
    isPatch = rhs.isPatch;
}
Sphere::~Sphere(){}
Sphere& Sphere::operator=(const Sphere& rhs){
    if(this != &rhs){
        Vector3d vector(rhs.coords);
        coords = vector;
        radius = rhs.radius;
        for (short i = 0; i < 8; i++)m_fill[i] = rhs.m_fill[i];
        isPatch = rhs.isPatch;
    }
    return *this;
}
Surface* Sphere::clone(){
    return new Sphere(*this);
}

//raySph
//checks to see if a sphere is hit
bool Sphere::intersect(const Ray& ray, double &min, double &max, Hit &hr){
    //solving for determinant
    //a is direction dot direction
    Vector3d vector = ray.eye - coords;
    Vector3d direction = ray.dir;
    double a = direction.squaredNorm();
    //b is direction dot eye - coords
    double b = direction.dot(vector);
    //c is eye - coords squared, minus radius squared
    double c = vector.squaredNorm()- (radius * radius);
    //first need to solve for discriminant to see if valid
    double discriminant = (b * b) - (a * c);

    //if discriminant is zero, then return false, else a sphere is hit
    if (discriminant < 0) return false;
    //else solve using quadratic and return true
    //changes the min and max
    double root1, root2;
    //roots calculated
    double ans = std::sqrt(discriminant);
    root1 = ((-1 * b) + ans)/(a);
    root2 = ((-1 * b) - ans)/(a);
    //choose appopriate root
    double t = ((root1 < 0) || ((root2 > 0) && (root2 < root1))) ? root2 : root1;
    //check that t is in the valid range
    if((t < min) || (t > max)) return false;
    max = t;
    hr.t = t;
    hr.inter = ray.eye + (t*ray.dir);
    hr.norm = ((hr.inter - coords) / radius);
    return true;
}

//getBounds
//based on radius of the sphere
void Sphere::getBounds(Vector3d& min, Vector3d& max) const{
    Vector3d r(radius, radius, radius);
    min = coords - r;
    max = coords + r;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// NODE CLASS ///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//Node constructor
//builds the node object
Node::Node(std::vector<Surface*>& objects, unsigned int start, unsigned int end){
    //initalize to prevent issues
    m_left = m_right = nullptr;
    m_min = m_max = Vector3d(0,0,0);
    m_data = nullptr; 

    //if count is equal to one, means base case and is the leaf node
    unsigned int count = end - start;
    if (count == 1) {
        m_data = objects[start];
        m_data->getBounds(m_min, m_max);
    }
    else{

        //initialize two vectors
        Vector3d cMin(1e9, 1e9, 1e9), cMax(-1e9, -1e9, -1e9);
        //then determine based on the box, what the min and max are of the bounding box
        for (unsigned int i = start; i < end; ++i) {
            Vector3d bMin, bMax;
            //gets bounds based on the object
            objects[i]->getBounds(bMin, bMax);
            Vector3d centroid = 0.5 * (bMin + bMax);
            //determines which to use
            for (short j = 0; j < 3; ++j) {
                cMin[j] = (cMin[j] < centroid[j]) ? cMin[j] : centroid[j];
                cMax[j] = (cMax[j] > centroid[j]) ? cMax[j] : centroid[j];
            }
        }
        //determines which axis to use based on which has longer 
        //in edge case where it need to go around, it is accounted for
        Vector3d spread = cMax - cMin;
        short axis = (spread[0] > spread[1] && spread[0] > spread[2]) ? 0 : ((spread[1] > spread[2]) ? 1 : 2);
        if (spread[axis] == 0) axis = (axis + 1) % 3; // handle flat axes

        //use quick select to get a rough estimate of a middle partition
        //done to attempt to save some runtime rather than sorting each time
        unsigned int mid = start + count / 2;
        quickselect(objects, start, end - 1, mid, axis);
        //now recursively create the nodes, then create the surrounding box
        m_left = new Node(objects, start, mid);
        m_right = new Node(objects, mid, end);
        surroundingBox(m_left->m_min, m_left->m_max, m_right->m_min, m_right->m_max, m_min, m_max);
    }

}

//boxIntersect
//check to see if ray intersects with bounding box
bool Node::boxIntersect(const Vector3d& min, const Vector3d& max, const Ray& ray, double tmin, double tmax){
    for (short i = 0; i < 3; i++) {
        //use the reciprocal of the direction and find intersection distances 
        double recip = 1.0 / ray.dir[i];
        double t0 = (min[i] - ray.eye[i]) * recip;
        double t1 = (max[i] - ray.eye[i]) * recip;
        if (recip < 0.0) std::swap(t0, t1);
        //find the lesser values
        tmin = (t0 > tmin) ? t0 : tmin;;
        tmax = (t1 > tmax) ? tmax : t1;

        //if the interval is logically inconsistent, return false
        if (tmax < tmin) return false;
    }
    //else if all 3 axes pass the test, return true
    return true;
}

//surroundingBox
//determines the bounding box
void Node::surroundingBox(const Vector3d& min1, const Vector3d& max1, const Vector3d& min2, const Vector3d& max2, Vector3d& outMin, Vector3d& outMax){
    //fmin and fmax used to prevent NaN errors
    for (short i = 0; i < 3; i++) {
        outMin[i] = std::fmin(min1[i], min2[i]);
        outMax[i] = std::fmax(max1[i], max2[i]);
    }
}

//compare
//compares centroids to see which is lesser
bool Node::compare(Surface* a, Surface *b, short axis){
    //vectors initialized
    Vector3d aMin, aMax, bMin, bMax;
    a->getBounds(aMin, aMax);
    b->getBounds(bMin, bMax);
    double aCenter = 0.5 * (aMin[axis] + aMax[axis]);
    double bCenter = 0.5 * (bMin[axis] + bMax[axis]);
    return aCenter < bCenter;
}

//nelement
//applies nth element algorithm
unsigned int Node::nelement(std::vector<Surface*>& surfaces, unsigned int start, unsigned int end, short axis){
    Surface* pivot = surfaces[end];
    int index = start;
    for(unsigned int i = start; i < end; ++i){
        if (compare(surfaces[i], pivot, axis)){
            std::swap(surfaces[i], surfaces[index]);
            ++index;
        }
    }
    //moves pivot to correct position
    std::swap(surfaces[index], surfaces[end]);
    return index;
}

//quick select
//attempts to use a simplified version of quicksort
void Node::quickselect(std::vector<Surface*>& surfaces, unsigned int start, unsigned int end, unsigned int mid, short axis){
    while (start <= end) {
        unsigned int pivot = nelement(surfaces, start, end, axis);
        if (pivot == mid) {
            return;
        } else if (pivot < mid) {
            start = pivot + 1;
        } else {
            end = pivot - 1;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// VIEWPOINT CLASS ///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//All of the setters for viewpoint
Viewpoint::Viewpoint(){}
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

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// LIGHT CLASS ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//creates a light
Light::Light(){}
Light::Light(double x, double y, double z){
    Vector3d vector(x, y, z);
    coords = vector; 
}