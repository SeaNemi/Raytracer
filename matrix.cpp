#include "matrix.h"

////////////////////////////////////////////////////////////
/////////////////// VECTOR3d CLASS ////////////////////////
//////////////////////////////////////////////////////////

//constructors
Vector3d::Vector3d(){
    m_x = m_y = m_z = 0;
}
Vector3d::Vector3d(double a, double b, double c){
    m_x = a;
    m_y = b;
    m_z = c;
}
Vector3d::Vector3d(const Vector3d& vec){
    m_x = vec.x();
    m_y = vec.y();
    m_z = vec.z();
}

//empty destructor since nothing is dynamically allocated
Vector3d::~Vector3d(){}

//cross product
Vector3d Vector3d::cross(const Vector3d& vec) const{
    return Vector3d((m_y*vec.z() - m_z*vec.y()),((m_z*vec.x() - m_x*vec.z())),(m_x*vec.y() - m_y*vec.x()));
}

//norm magntitude calculations
double Vector3d::norm() const{
    return std::sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
}
double Vector3d::magnitude() const{
    return norm();
}
double Vector3d::squaredNorm() const{
    return (m_x*m_x + m_y*m_y + m_z*m_z);
}
double Vector3d::dot(const Vector3d& vec) const{
    return ((m_x*vec.x())+ (m_y*vec.y()) + (m_z*vec.z()));
}

//normalize
//normalizes the vector
void Vector3d::normalize(){
    double nor = norm();
    if (nor == 0) return; //return to prevent divide by 0 error
    m_x /= nor;
    m_y /= nor;
    m_z /= nor;
}
//normalized
//returns a normalized copy of the vector
Vector3d Vector3d::normalized() const{
    double nor = norm();
    if (nor==0) return Vector3d();
    return Vector3d((m_x/nor), (m_y/nor), (m_z/nor));
}


//at and overloaded bracket operators
//return m_x m_y or m_z based on index, out of bounds not included as not necessary in this context
double& Vector3d::at(short index){
    if(index == 0) return m_x;
    else if(index == 1) return m_y;
    else return m_z;
}

double& Vector3d::operator[](short index){
    return at(index);
}
const double& Vector3d::operator[](short index) const{
    if(index == 0) return m_x;
    else if(index == 1) return m_y;
    else return m_z;
}

//overloaded operators
Vector3d Vector3d::operator+(const Vector3d& vec) const{
    return Vector3d ((m_x+vec.x()), (m_y+vec.y()), (m_z+vec.z()));
}

Vector3d Vector3d::operator-(const Vector3d& vec) const{
    return Vector3d ((m_x-vec.x()), (m_y-vec.y()), (m_z-vec.z()));
}
Vector3d Vector3d::operator*(double a) const{
    return Vector3d((m_x*a), (m_y*a), (m_z*a));
}
Vector3d Vector3d::operator/(double a) const{
    return Vector3d ((m_x/a), (m_y/a), (m_z/a));
}
Vector3d& Vector3d::operator+=(const Vector3d& vec){
    m_x += vec.x();
    m_y += vec.y();
    m_z += vec.z();
    return *this;
}
Vector3d& Vector3d::operator-=(const Vector3d& vec){
    m_x -= vec.x();
    m_y -= vec.y();
    m_z -= vec.z();
    return *this;   
}
Vector3d& Vector3d::operator*=(double a){
    m_x *= a;
    m_y *= a;
    m_z *= a;
    return *this;
}
Vector3d& Vector3d::operator/=(double a){
    m_x /= a;
    m_y /= a;
    m_z /= a;
    return *this;    
}

//assignment operator
Vector3d& Vector3d::operator=(const Vector3d& vec){
    if(this != &vec){
        m_x = vec.x();
        m_y = vec.y();
        m_z = vec.z();
    }
    return *this;
}


//boolean check
bool Vector3d::operator==(const Vector3d& vec) const{
    return ((m_x == vec.x()) && (m_y == vec.y()) && (m_z == vec.z()));
}
bool Vector3d::operator!=(const Vector3d& vec) const{
    return ((m_x != vec.x()) || (m_y != vec.y()) || (m_z != vec.z()));
}

//return member variables
double Vector3d::x() const{
    return m_x;
}
double Vector3d::y() const{
    return m_y;
}
double Vector3d::z() const{
    return m_z;
}

////////////////////////////////////////////////////////////
/////////////////// VECTOR2d CLASS ////////////////////////
//////////////////////////////////////////////////////////
//constructors
Vector2d::Vector2d(){
    m_x = m_y = 0;
}
Vector2d::Vector2d(double a, double b){
    m_x = a;
    m_y = b;
}
Vector2d::Vector2d(const Vector2d& vec){
    m_x = vec.x();
    m_y = vec.y();
}

//cross product
double Vector2d::cross(const Vector2d& vec) const{
    return (m_x*vec.y()- m_y*vec.x());
}

//dot product
double Vector2d::dot(const Vector2d& vec) const{
    return ((m_x*vec.x())+ (m_y*vec.y()));
}

//norm magntitude calculations
double Vector2d::norm() const{
    return std::sqrt(m_x*m_x + m_y*m_y);
}
double Vector2d::magnitude() const{
    return norm();
}
double Vector2d::squaredNorm() const{
    return (m_x*m_x + m_y*m_y);
}

Vector2d Vector2d::normalized() const{
    double nor = norm();
    if (nor==0) return Vector2d();
    return Vector2d((m_x/nor), (m_y/nor));
}

void Vector2d::normalize(){
    double nor = norm();
    if (nor == 0) return; //return to prevent divide by 0 error
    m_x /= nor;
    m_y /= nor;    
}


//empty destructor since nothing is dynamically allocated
Vector2d::~Vector2d(){}

double& Vector2d::at(short index){
    if(index == 0) return m_x;
    else return m_y;
}

double& Vector2d::operator[](short index){
    return at(index);
}
const double& Vector2d::operator[](short index) const{
    if(index == 0) return m_x;
    else return m_y;
}

//overloaded operators
Vector2d Vector2d::operator+(const Vector2d& vec) const{
    return Vector2d((m_x+vec.x()), (m_y+vec.y()));
}

Vector2d Vector2d::operator-(const Vector2d& vec) const{
    return Vector2d((m_x-vec.x()), (m_y-vec.y()));
}
Vector2d Vector2d::operator*(double a) const{
    return Vector2d((m_x*a), (m_y*a));
}
Vector2d Vector2d::operator/(double a) const{
    return Vector2d((m_x/a), (m_y/a));
}
Vector2d& Vector2d::operator+=(const Vector2d& vec){
    m_x += vec.x();
    m_y += vec.y();
    return *this;
}
Vector2d& Vector2d::operator-=(const Vector2d& vec){
    m_x -= vec.x();
    m_y -= vec.y();
    return *this;   
}
Vector2d& Vector2d::operator*=(double a){
    m_x *= a;
    m_y *= a;
    return *this;
}
Vector2d& Vector2d::operator/=(double a){
    m_x /= a;
    m_y /= a;
    return *this;    
}

//assignment operator
Vector2d& Vector2d::operator=(const Vector2d& vec){
    if(this != &vec){
        m_x = vec.x();
        m_y = vec.y();
    }
    return *this;
}

//boolean check
bool Vector2d::operator==(const Vector2d& vec) const{
    return ((m_x == vec.x()) && (m_y == vec.y()));
}
bool Vector2d::operator!=(const Vector2d& vec) const{
    return ((m_x != vec.x()) || (m_y != vec.y()));
}

double Vector2d::x() const{
    return m_x;
}
double Vector2d::y() const{
    return m_y;
}

//overloaded constructor
Vector2d::Vector2d(const Vector3d& vec, short dir)
    : Vector2d(project(vec, dir)) {}


//project
//used to aid in projection of object
Vector2d Vector2d::project(const Vector3d &vec, short dir) {
    switch (dir) {
        case 0:
            return Vector2d(vec.y(),vec.z());
        case 1:
            return Vector2d(vec.x(),vec.z());
        case 2:
            return Vector2d(vec.x(),vec.y());
    }
    return Vector2d(1.0, 1.0);
}



////////////////////////////////////////////////////////////
/////////////////// Matrix3d FUNCTIONS ////////////////////////
//////////////////////////////////////////////////////////

//constructors
Matrix3d::Matrix3d(){
    for(short i = 0; i < 3; i++) for(short j = 0; j < 3; j++) m_mat[i][j] = 0;
}

Matrix3d::Matrix3d(double arr[3][3]){
    for(short i = 0; i < 3; i++) for(short j = 0; j < 3; j++) m_mat[i][j] = arr[i][j];
}

Matrix3d::Matrix3d(const Matrix3d& arr){
    for(short i = 0; i < 3; i++) for(short j = 0; j < 3; j++) m_mat[i][j] = arr(i, j);
}

Matrix3d::Matrix3d(double arr[9]){
    short t = 0;
    for(short i = 0; i < 3; i++){
        for(short j = 0; j < 3; j++){
            m_mat[i][j] = arr[t];
            t++;
        }
    }
}

Matrix3d::Matrix3d(double a, double b, double c, double d, double e, double f, double g, double h, double i){
    m_mat[0][0]= a;
    m_mat[0][1]= b;
    m_mat[0][2]= c;
    m_mat[1][0]= d;
    m_mat[1][1]= e;
    m_mat[1][2]= f;
    m_mat[2][0]= g;
    m_mat[2][1]= h;
    m_mat[2][2]= i;
}

Matrix3d::Matrix3d(const Vector3d& vay, const Vector3d& vee, const Vector3d& viu){
    m_mat[0][0]= vay.x();
    m_mat[0][1]= vee.x();
    m_mat[0][2]= viu.x();
    m_mat[1][0]= vay.y();
    m_mat[1][1]= vee.y();
    m_mat[1][2]= viu.y();
    m_mat[2][0]= vay.z();
    m_mat[2][1]= vee.z();
    m_mat[2][2]= viu.z();
}

//destructor, empty but clear function is included
Matrix3d::~Matrix3d(){}
void Matrix3d::clear(){
    for(short i = 0; i < 3; i++) for(short j = 0; j < 3; j++) m_mat[i][j] = 0;    
}

//indexing elements
double& Matrix3d::operator()(short i, short j){
    return m_mat[i][j];
}
double Matrix3d::operator()(short i, short j) const{
    return m_mat[i][j];
}
double& Matrix3d::at(short i, short j){
    return m_mat[i][j];
}

//Overloaded Assignment operator
Matrix3d& Matrix3d::operator=(const Matrix3d& arr){
    if(this != &arr){
        clear();
        for(short i = 0; i < 3; i++) for(short j = 0; j < 3; j++) m_mat[i][j] = arr(i, j);
    }
    return *this;
}


//Overloaded mathematics operator
Matrix3d Matrix3d::operator+(const Matrix3d& arr) const{
    Matrix3d ret;
    for(short i=0; i < 3; i++){
        for(short j =0; j < 3; j++){
            ret(i, j) = m_mat[i][j] + arr(i,j);
        }
    }
    return ret;
}
Matrix3d Matrix3d::operator-(const Matrix3d& arr) const{
    Matrix3d ret;
    for(short i=0; i < 3; i++){
        for(short j =0; j < 3; j++){
            ret(i, j) = m_mat[i][j] - arr(i,j);
        }
    }
    return ret;
}
Matrix3d Matrix3d::operator*(const Matrix3d& arr) const{
    Matrix3d ret;
    for (short i = 0; i < 3; i++) {
        for (short j = 0; j < 3; j++) {
            ret(i, j) = 0;
            for (short k = 0; k < 3; k++) {
                ret(i, j) += m_mat[i][k] * arr(k, j);
            }
        }
    }

    return ret;
}
Matrix3d Matrix3d::operator*(double a) const{
    Matrix3d ret;
    for(short i=0; i < 3; i++){
        for(short j =0; j < 3; j++){
            ret(i, j) = m_mat[i][j] * a;
        }
    }
    return ret;
}
Matrix3d Matrix3d::operator/(double a) const{
    Matrix3d ret;
    for(short i=0; i < 3; i++){
        for(short j =0; j < 3; j++){
            ret(i, j) = m_mat[i][j] / a;
        }
    }
    return ret;
}
Matrix3d& Matrix3d::operator*=(double a){
    for(short i = 0; i < 3; i++) for(short j = 0; j < 3; j++) m_mat[i][j] *= a;
    return *this;
}
Matrix3d& Matrix3d::operator/=(double a){
    for(short i = 0; i < 3; i++) for(short j = 0; j < 3; j++) m_mat[i][j] /= a;
    return *this;
}
Matrix3d& Matrix3d::operator+=(const Matrix3d& arr){
    for(short i = 0; i < 3; i++) for(short j = 0; j < 3; j++) m_mat[i][j] += arr(i, j);
    return *this;
}
Matrix3d& Matrix3d::operator-=(const Matrix3d& arr){
    for(short i = 0; i < 3; i++) for(short j = 0; j < 3; j++) m_mat[i][j] -= arr(i, j);
    return *this;
}

Vector3d Matrix3d::operator*(const Vector3d& vec) const{
    double x = m_mat[0][0] * vec.x() + m_mat[0][1] * vec.y() + m_mat[0][2] * vec.z();
    double y = m_mat[1][0] * vec.x() + m_mat[1][1] * vec.y() + m_mat[1][2] * vec.z();
    double z = m_mat[2][0] * vec.x() + m_mat[2][1] * vec.y() + m_mat[2][2] * vec.z();
    return Vector3d(x, y, z);
}

//transpose
Matrix3d Matrix3d::transpose() const{
    Matrix3d copy;
    for(short i = 0; i < 3; i++) for(short j = 0; j < 3; j++) copy(j, i) = m_mat[i][j];
    return copy;
}

double Matrix3d::det() const{
    return determinant();
}

double Matrix3d::determinant() const{
    double start = m_mat[0][0] * (m_mat[1][1]*m_mat[2][2] - m_mat[1][2]*m_mat[2][1]);
    start -= m_mat[0][1] * (m_mat[1][0]*m_mat[2][2] - m_mat[2][0]*m_mat[1][2]);
    start += m_mat[0][2] * (m_mat[1][0]*m_mat[2][1] - m_mat[2][0]*m_mat[1][1]);
    return start;
}

//row and col
void Matrix3d::row(short index, const Vector3d& vec){
    if(index < 0 || index > 2) return;
    for(short i = 0; i < 3; i++) m_mat[index][i] = vec[i];
}
void Matrix3d::col(short index, const Vector3d& vec){
    if(index < 0 || index > 2) return;
    for(short i = 0; i < 3; i++) m_mat[i][index] = vec[i];
}

////////////////////////////////////////////////////////////
/////////////////// OUTSIDE FUNCTIONS ////////////////////////
//////////////////////////////////////////////////////////

Vector3d operator*(double a, const Vector3d& vec) {
    return vec * a;
}

Vector3d cross(const Vector3d& x, const Vector3d& y){
    return x.cross(y);
}