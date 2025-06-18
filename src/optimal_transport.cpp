#include "optimal_transport.h"
#include <iostream> 

static lbfgsfloatval_t evaluate(void* instance, const lbfgsfloatval_t* w, lbfgsfloatval_t* grad, const int n, const lbfgsfloatval_t step) {
    
    OptimalTransport* ot = static_cast<OptimalTransport*>(instance);
    for (int i = 0; i < ot->num_fluid_particles; i++) {
        ot->diagram.weights[i] = w[i];
    }
    ot->diagram.weights[ot->num_fluid_particles] = w[ot->num_fluid_particles];
    ot->diagram.compute();
    
    lbfgsfloatval_t res = 0.0;
    double fluid_vol = 0.0;
    
    for (int i = 0; i < ot->num_fluid_particles; i++) {
        double area = ot->diagram.cells[i].area();
        double desired_fluid_volume = ot->target_fluid_masses[i];
         // we compute first the integral ||x - yi||^2 over power cell
        res += ot->diagram.cells[i].compute_integral(ot->diagram.points[i]);
         // we subtract the Lagrange multiplier term
        res -= w[i] * (area - desired_fluid_volume);
        // gradient for fluid 
        grad[i] = -(desired_fluid_volume - area);
        fluid_vol += area;
    }
    
    double air_vol = 1.0 - fluid_vol;
    res -= w[ot->num_fluid_particles] * (air_vol - ot->desired_vol_air);
    // gradient for air 
    grad[ot->num_fluid_particles] = -(ot->desired_vol_air - air_vol);
    
    return -res;
}

//both progress and optimize were done with the helped of LLMs 
static int progress(void* instance, const lbfgsfloatval_t* x, const lbfgsfloatval_t* g, const lbfgsfloatval_t fx,
                    const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls) {
    return 0; 
}


void OptimalTransport::optimize() {
    int N_lbfgs_vars = num_fluid_particles + 1;
    lbfgsfloatval_t fx;
    std::vector<lbfgsfloatval_t> weights_lbfgs(N_lbfgs_vars);
    for (int i = 0; i < num_fluid_particles; i++) {
        weights_lbfgs[i] = diagram.weights[i];
    }
    weights_lbfgs[num_fluid_particles] = diagram.weights[num_fluid_particles];
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param); 
    int ret = lbfgs(N_lbfgs_vars, &weights_lbfgs[0], &fx, evaluate, progress, static_cast<void*>(this), &param);
    // std::cout << "L-BFGS optimization terminated with status code " << ret << std::endl;
    for (int i = 0; i < num_fluid_particles; i++) {
        diagram.weights[i] = weights_lbfgs[i];
    }
    diagram.weights[num_fluid_particles] = weights_lbfgs[num_fluid_particles];
    diagram.compute();
}