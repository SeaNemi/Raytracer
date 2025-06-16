
#include <cmath>

class Matrix
{


};


class Vector3d
{
    public:

    //constructors and destructor
    Vector3d();
    Vector3d(double, double, double);
    Vector3d(const Vector3d&);
    ~Vector3d();

    //cross
    Vector3d cross(const Vector3d&) const;

    //Calculations for norm/magnitude and dot
    double norm() const;
    double magnitude() const;
    double normSquared() const;
    double dot(const Vector3d&) const;

    //normalize and normalized
    void normalize(); //normalizes vector
    Vector3d normalized() const; //returns normalized copy

    //return indexing
    double& at(int);
    double& operator[](int);


    //overloaded operators
    Vector3d operator+(const Vector3d&) const;
    Vector3d operator-(const Vector3d&) const;
    Vector3d operator*(double) const;
    Vector3d operator/(double) const;
    Vector3d& operator+=(const Vector3d&);
    Vector3d& operator-=(const Vector3d&);
    Vector3d& operator*=(double);
    Vector3d& operator/=(double);
    Vector3d& operator=(const Vector3d&);
    bool operator==(const Vector3d&);


    //member variables
    double x;
    double y;
    double z;
};