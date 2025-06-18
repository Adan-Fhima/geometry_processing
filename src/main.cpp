#include <cstdlib> 
#include "simulation.h" 

int main() {
    system("rm -f fluid_sim_*.png");
    simulate_fluid();
    return 0;
}