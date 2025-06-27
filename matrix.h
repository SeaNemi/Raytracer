//William Whatley
//matrix header file

//Header guards
#ifndef MATRIX_H
#define MATRIX_H

//library included and identity matrix
#include <cmath>
const double IDENTITY[3][3]={
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1}
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
    double squaredNorm() const;
    double dot(const Vector3d&) const;

    //normalize and normalized
    void normalize(); //normalizes vector
    Vector3d normalized() const; //returns normalized copy

    //return indexing
    double& at(short);
    double& operator[](short);
    const double& operator[](short) const;


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
    bool operator==(const Vector3d&) const;
    bool operator!=(const Vector3d&) const;

    //return member variables
    double x() const;
    double y() const;
    double z() const;

    //member variables
    private:
    double m_x;
    double m_y;
    double m_z;
};

class Vector2d
{
    public:
    //constructors and destructor
    Vector2d();
    Vector2d(double, double);
    Vector2d(const Vector2d&);
    ~Vector2d();

    //cross product
    double cross(const Vector2d& vec) const;

    //mathematical operations
    double dot(const Vector2d&) const;
    double norm() const;
    double magnitude() const;
    double squaredNorm() const;
    Vector2d normalized() const;
    void normalize();

    //return indexing
    double& at(short);
    double& operator[](short);
    const double& operator[](short) const;

    //overloaded operators
    Vector2d operator+(const Vector2d&) const;
    Vector2d operator-(const Vector2d&) const;
    Vector2d operator*(double) const;
    Vector2d operator/(double) const;
    Vector2d& operator+=(const Vector2d&);
    Vector2d& operator-=(const Vector2d&);
    Vector2d& operator*=(double);
    Vector2d& operator/=(double);
    Vector2d& operator=(const Vector2d&);
    bool operator==(const Vector2d&) const;
    bool operator!=(const Vector2d&) const;

    //return member variables
    double x() const;
    double y() const;

    //computer graphics specific code
    //overloaded constructor based on projection
    Vector2d(const Vector3d&, short);
    Vector2d project(const Vector3d&, short);

    //member variables
    private:
    double m_x;
    double m_y;
};

class Matrix3d
{
    public:
    Matrix3d();
    Matrix3d(double[3][3]);
    Matrix3d(double[9]);
    Matrix3d(double, double, double, double, double, double, double, double, double);
    Matrix3d(const Matrix3d&);
    Matrix3d(const Vector3d&, const Vector3d&, const Vector3d&);
    ~Matrix3d();
    void clear();

    //return indexing
    double& operator()(short, short);
    double operator()(short, short) const;
    double& at(short i, short j);

    //assignment operator
    Matrix3d& operator=(const Matrix3d&);

    //overloaded math constructors
    Matrix3d operator+(const Matrix3d&) const;
    Matrix3d operator-(const Matrix3d&) const;
    Matrix3d operator*(const Matrix3d&) const;
    Matrix3d operator*(double) const;    
    Matrix3d operator/(double) const;
    Matrix3d& operator*=(double);
    Matrix3d& operator/=(double);   
    Matrix3d& operator+=(const Matrix3d&);
    Matrix3d& operator-=(const Matrix3d&);

    //vector multiplication operator
    Vector3d operator*(const Vector3d&) const;

    //normal functions
    Matrix3d transpose() const;
    double det() const;
    double determinant() const;
    Matrix3d inverse() const;

    void row(short, const Vector3d&);
    void col(short, const Vector3d&);

    private:
    double m_mat[3][3];
};

//overloaded operator to deal with inconsistency issues with syntax
Vector3d operator*(double a, const Vector3d& v);
Vector3d cross(const Vector3d& x, const Vector3d& y);

#endif