#include "matrix.h"

//constructors
Vector3d::Vector3d(){
    x = y = z = 0;
}
Vector3d::Vector3d(double a, double b, double c){
    x = a;
    y = b;
    z = c;
}
Vector3d::Vector3d(const Vector3d& vec){
    x = vec.x;
    y = vec.y;
    z = vec.z;
}

//empty destructor since nvecing dynamically allocated
Vector3d::~Vector3d(){}

//cross product
Vector3d Vector3d::cross(const Vector3d& vec) const{
    return Vector3d((y*vec.z - z*vec.y),(-1*(z*vec.x - x*vec.z)),(x*vec.y-y*vec.x));
}

//norm magntitude calculations
double Vector3d::norm() const{
    return std::sqrt(x*x + y*y + z*z);
}
double Vector3d::magnitude() const{
    return norm();
}
double Vector3d::normSquared() const{
    return (x*x + y*y + z*z);
}
double Vector3d::dot(const Vector3d& vec) const{
    return ((x*vec.x)+ (y*vec.y) + (z*vec.z));
}

//normalize
//normalizes the vector
void Vector3d::normalize(){
    double nor = norm();
    if (nor == 0) return; //return to prevent divide by 0 error
    x /= nor;
    y /= nor;
    z /= nor;
}
//normalized
//returns a normalized copy of the vector
Vector3d Vector3d::normalized() const{
    double nor = norm();
    if (nor==0) return Vector3d();
    return Vector3d((x/nor), (y/nor), (z/nor));
}


//at and overloaded bracket operators
//return x y or z based on index, out of bounds not included as not necessary in this context
double& Vector3d::at(int index){
    if(index == 0) return x;
    else if(index == 1) return y;
    else return z;
}

double& Vector3d::operator[](int index){
    return at(index);
}

//overloaded operators
Vector3d Vector3d::operator+(const Vector3d& vec) const{
    return Vector3d ((x+vec.x), (y+vec.y), (z+vec.z));
}

Vector3d Vector3d::operator-(const Vector3d& vec) const{
    return Vector3d ((x-vec.x), (y-vec.y), (z-vec.z));
}
Vector3d Vector3d::operator*(double a) const{
    return Vector3d ((x*a), (y*a), (z*a));
}
Vector3d Vector3d::operator/(double a) const{
    return Vector3d ((x*a), (y*a), (z*a));
}
Vector3d& Vector3d::operator+=(const Vector3d& vec){
    x += vec.x;
    y += vec.y;
    z += vec.z;
    return *this;
}
Vector3d& Vector3d::operator-=(const Vector3d& vec){
    x -= vec.x;
    y -= vec.y;
    z -= vec.z;
    return *this;   
}
Vector3d& Vector3d::operator*=(double a){
    x *= a;
    y *= a;
    z *= a;
    return *this;
}
Vector3d& Vector3d::operator/=(double a){
    x /= a;
    y /= a;
    z /= a;
    return *this;    
}

Vector3d& Vector3d::operator=(const Vector3d& vec){
    if(this != &vec){
        x = vec.x;
        y = vec.y;
        z = vec.z;
    }
    return *this;
}      