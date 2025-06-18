#ifndef UTILS_H
#define UTILS_H

#include <random>
#include <omp.h> 
#include <cmath> 
#include <vector>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


extern thread_local std::default_random_engine engine;
extern thread_local std::uniform_real_distribution<double> uniform;


int sgn(double val);
class Vector {
public:
    double data[3];
    explicit Vector(double x = 0, double y = 0, double z = 0);
    double norm2() const;
    double norm() const;
    void normalize();
    double operator[](int i) const;
    double& operator[](int i);
};

Vector operator+(const Vector& a, const Vector& b);
Vector operator-(const Vector& a, const Vector& b);
Vector operator-(const Vector& a);
Vector operator*(double a, const Vector& b);
Vector operator*(const Vector& a, const Vector& b);
Vector operator*(const Vector& a, double b);
Vector operator/(const Vector& a, double b);
double dot(const Vector& a, const Vector& b);
Vector cross(const Vector& a, const Vector& b);

#endif // UTILS_H