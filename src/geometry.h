#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include "utils.h" 


class Polygon {
public:
    std::vector<Vector> vertices;

    double area() const;
    Vector centroid();
    double compute_integral(const Vector& Pi) const;
};


enum ParticleType {FLUID, AIR};

struct Particle {
    ParticleType type;
    Vector position;
    Vector velocity;
    double mass;
    double weight;
    Vector last_centroid;
};


class PowerDiagram {
public:
    PowerDiagram() = default;
    std::vector<Vector> points;
    std::vector<double> weights;
    std::vector<ParticleType> types;
    int num_fluid = 0;
    std::vector<Polygon> cells; 
    void compute(); 

private:
    bool isInside(const Vector& X, const Vector& P0, double w0, const Vector& Pi, double wi);
    Polygon clip_by_bisector(const Polygon& init_pol, const Vector& P0, double w0, const Vector& Pi, double wi);
    Polygon clip_by_disk(const Polygon& init_pol, const Vector& center, double radius_sq);
};

#endif // GEOMETRY_H