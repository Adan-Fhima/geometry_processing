#ifndef OPTIMAL_TRANSPORT_H
#define OPTIMAL_TRANSPORT_H

#include <vector>
#include "geometry.h"
#include "../lbfgs.h"    


static lbfgsfloatval_t evaluate(void* instance, const lbfgsfloatval_t* x, lbfgsfloatval_t* g, const int n, const lbfgsfloatval_t step);
static int progress(void* instance, const lbfgsfloatval_t* x, const lbfgsfloatval_t* g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls);


class OptimalTransport {
public:
    OptimalTransport() = default;
    void optimize();
    PowerDiagram diagram;
    std::vector<double> target_fluid_masses; 
    double desired_vol_air = 0.0; 
    int num_fluid_particles = 0;
    int num_air_particles = 0;
};

#endif // OPTIMAL_TRANSPORT_H