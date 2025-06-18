#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <string>

#include "geometry.h"
#include "utils.h"   


struct SimulationParameters {
    static constexpr double GRAVITY_ACCEL = -9.81;
    static constexpr double EPSILON_SQ = 0.004 * 0.004;
    static constexpr double DT = 0.01;
    static constexpr double MASS = 200.0;
    static constexpr int N_FLUID = 50;
    static constexpr int N_AIR = 400;
    static constexpr int NUM_FRAMES = 100;
    static constexpr double FLUID_RADIUS = 0.2;
};

void save_svg(const std::vector<Polygon>& polygons, std::string filename, const std::vector<Vector>* points = nullptr);
void save_frame(const std::vector<Polygon>& cells, const std::vector<Particle>& particles, std::string filename, int frameid);
void simulate_fluid();

#endif // SIMULATION_H