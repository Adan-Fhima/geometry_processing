#include "geometry.h"
#include <algorithm> 
#include <cmath>     

double Polygon::area() const {
    if (vertices.size() < 3) return 0.0;
    double result = 0.0;
    for (size_t i = 0, j = vertices.size() - 1; i < vertices.size(); j = i++) {
        result += (vertices[j][0] * vertices[i][1] - vertices[j][1] * vertices[i][0]);
    }
    return std::abs(result) * 0.5;
}

Vector Polygon::centroid() {
    if (vertices.size() < 3) return Vector(0., 0., 0.);
    Vector cent(0., 0., 0.);
    for (size_t i = 0; i < vertices.size(); i++) {
        size_t ip = i + 1;
        if (ip == vertices.size()) { ip = 0; }
        double crossP = (vertices[i][0] * vertices[ip][1] - vertices[ip][0] * vertices[i][1]);
        cent = cent + (vertices[i] + vertices[ip]) * crossP;
    }
    double pol_area = area();
    if (std::abs(pol_area) < 1e-10) return Vector(0., 0., 0.); 
    cent = cent / (6. * pol_area);
    return cent;
}

double Polygon::compute_integral(const Vector& Pi) const {
    if (vertices.size() < 3) return 0.0;
    double result = 0.0;
    for (size_t t = 0; t < vertices.size() - 2; t++) {
        Vector v0 = vertices[0];
        Vector v1 = vertices[t+1];
        Vector v2 = vertices[t+2];
        Vector d0 = v0 - Pi;
        Vector d1 = v1 - Pi;
        Vector d2 = v2 - Pi;
        Vector d[3] = {d0, d1, d2};
        double sum_res = 0.0;
        for (int k = 0; k < 3; k++) {
            for (int l = k; l < 3; l++) {
                sum_res += dot(d[k], d[l]);
            }
        }
        Vector c2c1 = v1 - v0;
        Vector c3c1 = v2 - v0;
        double areaT = 0.5 * std::abs(c2c1[0] * c3c1[1] - c2c1[1] * c3c1[0]);
        result += (areaT / 6.0) * sum_res;
    }
    return result;
}


bool PowerDiagram::isInside(const Vector& X, const Vector& P0, double w0, const Vector& Pi, double wi) {
    return (X - P0).norm2() - w0 <= (X - Pi).norm2() - wi;
}

Polygon PowerDiagram::clip_by_bisector(const Polygon& init_pol, const Vector& P0, double w0, const Vector& Pi, double wi) {
    Polygon result;
    size_t n = init_pol.vertices.size();
    if (n == 0) return result;

    for (size_t i = 0; i < n; i++) {
        Vector A = init_pol.vertices[i];
        Vector B = init_pol.vertices[(i + 1) % n];
        bool A_inside = isInside(A, P0, w0, Pi, wi);
        bool B_inside = isInside(B, P0, w0, Pi, wi);
        Vector midPoint = 0.5 * (P0 + Pi);
        double normSquared = (P0 - Pi).norm2();
        Vector M_prime = midPoint + ((w0 - wi) / (2.0 * normSquared)) * (Pi - P0);

        if (B_inside) {
            if (!A_inside) {
                double denominator = dot(B - A, Pi - P0);
                if (std::abs(denominator) > 1e-10) {
                    double t = dot(M_prime - A, Pi - P0) / denominator;
                    Vector P = A + t * (B - A);
                    result.vertices.push_back(P);
                }
            }
            result.vertices.push_back(B);
        } else if (A_inside) {
            double denominator = dot(B - A, Pi - P0);
            if (std::abs(denominator) > 1e-10) {
                double t = dot(M_prime - A, Pi - P0) / denominator;
                Vector P = A + t * (B - A);
                result.vertices.push_back(P);
            }
        }
    }
    return result;
}

Polygon PowerDiagram::clip_by_disk(const Polygon& init_pol, const Vector& center, double radius_sq) {
    Polygon result;
    size_t n = init_pol.vertices.size();
    if (n == 0) return result;

    for (size_t i = 0; i < n; i++) {
        Vector A = init_pol.vertices[i];
        Vector B = init_pol.vertices[(i + 1) % n];
        bool A_inside = (A - center).norm2() <= radius_sq;
        bool B_inside = (B - center).norm2() <= radius_sq;

        if (B_inside) {
            if (!A_inside) {
                Vector dp = B - A;
                Vector ac = A - center;
                double a = dp.norm2();
                double b = 2.0 * dot(ac, dp);
                double c = ac.norm2() - radius_sq;
                
                double discriminant = b * b - 4 * a * c;
                if (discriminant >= 0) {
                    double t = (-b - std::sqrt(discriminant)) / (2.0 * a);
                    if (t >= 0.0 && t <= 1.0) {
                        result.vertices.push_back(A + t * dp);
                    }
                }
            }
            result.vertices.push_back(B);
        } else if (A_inside) {
            Vector dp = B - A;
            Vector ac = A - center;
            double a = dp.norm2();
            double b = 2.0 * dot(ac, dp);
            double c = ac.norm2() - radius_sq;
            
            double discr = b * b - 4 * a * c;
            if (discr >= 0) {
                double t = (-b + std::sqrt(discr)) / (2.0 * a);
                if (t >= 0.0 && t <= 1.0) {
                    result.vertices.push_back(A + t * dp);
                }
            }
        }
    }
    return result;
}
void PowerDiagram::compute() {
    Polygon Square; 
    Square.vertices.push_back(Vector(0, 0, 0));
    Square.vertices.push_back(Vector(1, 0, 0));
    Square.vertices.push_back(Vector(1, 1, 0));
    Square.vertices.push_back(Vector(0, 1, 0));

    cells.resize(points.size()); 

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < points.size(); i++) {
        Polygon V = Square;
        Vector Pi = points[i];
        double wi = weights[i];
        for (size_t j = 0; j < points.size(); j++) {
            if (i == j) continue; 
            if (V.vertices.empty()) break;
        
            Vector Pj = points[j];
            double wj = weights[j];
            V = clip_by_bisector(V, Pi, wi, Pj, wj);
        }

        if (types[i] == FLUID) {
            double w_air = weights[num_fluid];
            double radius_sq = wi - w_air;
            if (radius_sq < 0) radius_sq = 0;

            if (!V.vertices.empty()) {
                V = clip_by_disk(V, Pi, radius_sq);
            }
        }
        cells[i] = V;
    }
}