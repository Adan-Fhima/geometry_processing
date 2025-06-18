#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <omp.h>
#include <algorithm>
#include <cstring>
#include "lbfgs.h"

thread_local std::default_random_engine engine(omp_get_thread_num() + 1); 
thread_local std::uniform_real_distribution<double> uniform(0.0, 1.0);

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        if (n > 1e-10) {
            data[0] /= n;
            data[1] /= n;
            data[2] /= n;
        }
    }
    double operator[](int i) const { return data[i]; };
    double &operator[](int i) { return data[i]; };
    double data[3];
};

Vector operator+(const Vector &a, const Vector &b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector &a, const Vector &b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator-(const Vector &a) {
    return Vector(-a[0], -a[1], -a[2]);
}

Vector operator*(const double a, const Vector &b) {
    return Vector(a * b[0], a * b[1], a * b[2]);
}

Vector operator*(const Vector &a, const Vector &b) {
  return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}

Vector operator*(const Vector &a, const double b) {
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}
Vector operator/(const Vector &a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector &a, const Vector &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector &a, const Vector &b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Polygon {  
public:
    std::vector<Vector> vertices;
    
    double area() const {
        if (vertices.size() < 3) return 0.0;
        double result = 0.0;
        for (int i = 0, j = vertices.size() - 1; i < vertices.size(); j = i++) {
            result += (vertices[j][0] * vertices[i][1] - vertices[j][1] * vertices[i][0]);
        }
        return std::abs(result) * 0.5;
    }
    
    Vector centroid(){
        if (vertices.size() <3) return Vector(0.,0.,0.);
        Vector cent(0.,0.,0.);
        for (int i = 0; i< vertices.size(); i++){
            int ip = i + 1;
            if (ip == vertices.size()) {ip = 0;}
            double crossP = (vertices[i][0] * vertices[ip][1] - vertices[ip][0] * vertices[i][1]);
            cent = cent + (vertices[i] + vertices[ip]) * crossP;
        }
        double pol_area = area();
        if (std::abs(pol_area)<1e-10) return Vector(0.,0.,0.);
        cent = cent / (6.* pol_area);
        return cent;
    }
    //compute int (x - pi).norm2 
    double compute_integral(const Vector& Pi) const {   
        if (vertices.size() < 3) return 0.0;
        double result = 0.0;
        for (int t = 0; t < vertices.size() - 2; t++) {
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

};



class PowerDiagram {
public:
    PowerDiagram() {};

     bool isInside(const Vector& X, const Vector& P0, double w0, const Vector& Pi, double wi) {
        return (X - P0).norm2() - w0 <= (X - Pi).norm2() - wi;
    }

   // Sutherland-Hodgman algorithm

    Polygon clip_by_bisector(const Polygon& init_pol, const Vector& P0, double w0, const Vector& Pi, double wi) {
        Polygon result;
        int n = init_pol.vertices.size();
        for (int i = 0; i < n; i++) {
            Vector A = init_pol.vertices[(i == 0) ? n - 1 : i - 1]; 
            Vector B = init_pol.vertices[i];                         
            bool A_inside = isInside(A, P0, w0, Pi, wi);
            bool B_inside = isInside(B, P0, w0, Pi, wi);
            Vector midPoint = 0.5 * (P0 + Pi);
            double normSquared = (P0 - Pi).norm2();
            Vector M_prime = midPoint + ((w0 - wi) / (2.0 * normSquared)) * (Pi - P0);
            
            if (B_inside) {
                if (!A_inside) {
                    double t = dot(M_prime - A, Pi - P0) / dot(B - A, Pi - P0);
                    Vector P = A + t * (B - A);
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            }
            else if (A_inside) {
                double t = dot(M_prime - A, Pi - P0) / dot(B - A, Pi - P0);
                Vector P = A + t * (B - A);
                result.vertices.push_back(P);
            }
        }
        return result;
    }


    void compute() {
        Polygon Square;
        Square.vertices.push_back(Vector(0, 0, 0));
        Square.vertices.push_back(Vector(1, 0, 0));
        Square.vertices.push_back(Vector(1, 1, 0));
        Square.vertices.push_back(Vector(0, 1, 0));
        
        cells.resize(points.size());
        
        for (size_t i = 0; i < points.size(); i++) {
            Polygon V = Square;
            for (size_t j = 0; j < points.size(); j++) {
                if (i == j) continue;
                if (V.vertices.empty()) break;
                V = clip_by_bisector(V, points[i], weights[i], points[j], weights[j]);
            }
            cells[i] = V;
        }
    }

    std::vector<Vector> points;
    std::vector<double> weights;
    std::vector<Polygon> cells;
};

class OptimalTransport {
public:
    OptimalTransport() {};
    void optimize();
    PowerDiagram diagram;
    std::vector<double> target_masses;
};

static lbfgsfloatval_t evaluate(void* instance, const lbfgsfloatval_t* x, lbfgsfloatval_t* g, const int n, const lbfgsfloatval_t step) {
    OptimalTransport* ot = (OptimalTransport*)(instance);
    for (int i = 0; i < n; i++) {
        ot->diagram.weights[i] = x[i];
    }
    ot->diagram.compute();
    
    //Section
    lbfgsfloatval_t fx = 0.0;    
    for (int i = 0; i < n; i++) {
        double pol_area = ot->diagram.cells[i].area();
        double target_area = ot->target_masses[i];
        fx += ot->diagram.cells[i].compute_integral(ot->diagram.points[i]);
        fx -= x[i] * (pol_area - target_area);
        g[i] = -(target_area - pol_area);
    }
    return -fx; 
}

static int progress(void* instance, const lbfgsfloatval_t* x, const lbfgsfloatval_t* g, const lbfgsfloatval_t fx, 
                    const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls) {
    return 0;
}

void OptimalTransport::optimize() {
    int N = diagram.weights.size();
    lbfgsfloatval_t fx;
    std::vector<lbfgsfloatval_t> weights(N);
    for (int i = 0; i < N; i++) {
        weights[i] = diagram.weights[i];
    }

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    param.max_iterations = 1000;
    param.epsilon = 1e-6;
    int ret = lbfgs(N, &weights[0], &fx, evaluate, progress, (void*)this, &param);
    printf("L-BFGS optimization terminated with status code %d\n", ret);
    for (int i = 0; i < N; i++) {
        diagram.weights[i] = weights[i];
    }
    diagram.compute();
}

void save_svg(const std::vector<Polygon>& polygons, std::string filename, const std::vector<Vector>* points = NULL) {
    FILE* f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (size_t i = 0; i < polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (size_t j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"none\" stroke = \"black\" stroke-width=\"1\"/>\n");
        fprintf(f, "</g>\n");
    }
    if (points) {
        fprintf(f, "<g>\n");        
        for (size_t i = 0; i < points->size(); i++) {
            fprintf(f, "<circle cx = \"%3.3f\" cy = \"%3.3f\" r = \"2\" fill=\"black\"/>\n", 
                    (*points)[i][0]*1000., 1000.-(*points)[i][1]*1000);
        }
        fprintf(f, "</g>\n");
    }

    fprintf(f, "</svg>\n");
    fclose(f);
}

int main() {
    int N = 200;
    OptimalTransport ot;
    for (int i = 0; i < N; i++) {
        ot.diagram.points.push_back(Vector(uniform(engine), uniform(engine), 0.0));
        ot.diagram.weights.push_back(0.0); 
    }
    ot.target_masses.resize(N, 1.0 / N);

    ot.diagram.compute();
    save_svg(ot.diagram.cells, "voronoi_standard.svg", &ot.diagram.points);
    
    //OT Voronoi
    printf("Optimizing for uniform distribution...\n");
    ot.optimize();
    double total_area = 0.0;
    double max_area = 0.0, min_area = 1.0;
    for (int i = 0; i < N; i++) {
        double area = ot.diagram.cells[i].area();
        // printf("area cell is %.6f", area, "\n");
        total_area += area;
        max_area = std::max(max_area, area);
        min_area = std::min(min_area, area);
    }
    printf("Total area: %.6f\n", total_area);
    printf("Area range: [%.6f, %.6f] (target: %.6f)\n", min_area, max_area, 1.0/N);
    printf("Area variation: %.2f%%\n", 100.0 * (max_area - min_area) / (1.0/N));
    save_svg(ot.diagram.cells, "power_diagram_optimized.svg", &ot.diagram.points);
    
    return 0;
}