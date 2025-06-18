#include "simulation.h"
#include "optimal_transport.h" 
#include "utils.h"         
#include <iostream>          
#include <sstream>        
#include <iomanip>            


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb_image_write.h"


void save_svg(const std::vector<Polygon>& polygons, std::string filename, const std::vector<Vector>* points) {
    FILE* f = std::fopen(filename.c_str(), "w+");
    if (!f) {
        std::cerr << "Error: Could not open SVG file " << filename << std::endl;
        return;
    }
    std::fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (size_t i = 0; i < polygons.size(); i++) {
        std::fprintf(f, "<g>\n");
        std::fprintf(f, "<polygon points = \"");
        for (size_t j = 0; j < polygons[i].vertices.size(); j++) {
            std::fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        std::fprintf(f, "\"\nfill = \"none\" stroke = \"black\" stroke-width=\"1\"/>\n");
        std::fprintf(f, "</g>\n");
    }
    if (points) {
        std::fprintf(f, "<g>\n");
        for (size_t i = 0; i < points->size(); i++) {
            std::fprintf(f, "<circle cx = \"%3.3f\" cy = \"%3.3f\" r = \"2\" fill=\"black\"/>\n",
                    (*points)[i][0]*1000., 1000.-(*points)[i][1]*1000);
        }
        std::fprintf(f, "</g>\n");
    }
    std::fprintf(f, "</svg>\n");
    std::fclose(f);
}

void save_frame(const std::vector<Polygon>& cells, const std::vector<Particle>& particles, std::string filename, int frameid) {
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W * H * 3, 255); 

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < cells.size(); i++) {
        double bminx = 1e9, bminy = 1e9, bmaxx = -1e9, bmaxy = -1e9;
        for (size_t j = 0; j < cells[i].vertices.size(); j++) {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(static_cast<double>(W - 1), std::max(0., W * bminx));
        bminy = std::min(static_cast<double>(H - 1), std::max(0., H * bminy));
        bmaxx = std::min(static_cast<double>(W - 1), std::max(0., W * bmaxx));
        bmaxy = std::min(static_cast<double>(H - 1), std::max(0., H * bmaxy));

        int min_x_int = static_cast<int>(bminx);
        int max_x_int = static_cast<int>(bmaxx) + 1;
        int min_y_int = static_cast<int>(bminy);
        int max_y_int = static_cast<int>(bmaxy) + 1;

        for (int y = min_y_int; y < max_y_int; y++) {
            for (int x = min_x_int; x < max_x_int; x++) {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1e9; 
                for (size_t j = 0; j < cells[i].vertices.size(); j++) {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                    double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                    double det = (x - x0) * (y1 - y0) - (y - y0) * (x1 - x0); 

                    int sign = sgn(det);
                    if (prevSign == 0) prevSign = sign;
                    else if (sign == 0) sign = prevSign;
                    else if (sign != prevSign) { 
                        isInside = false;
                        break;
                    }
                    prevSign = sign;
                    double edgeLen_sq = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0);
                    double distEdge = 1e9;
                    if (edgeLen_sq > 1e-10) { 
                        distEdge = std::abs(det) / std::sqrt(edgeLen_sq);
                        double dotp_proj = dot(Vector(x - x0, y - y0, 0), Vector(x1 - x0, y1 - y0, 0));
                        if (dotp_proj < 0 || dotp_proj > edgeLen_sq) distEdge = 1e9;
                    }
                    mindistEdge = std::min(mindistEdge, distEdge);
                }

                if (isInside) {
                    if (static_cast<size_t>(i) < particles.size() && particles[i].type == FLUID) {
                        image[((H - y - 1) * W + x) * 3] = 0;
                        image[((H - y - 1) * W + x) * 3 + 1] = 0;
                        image[((H - y - 1) * W + x) * 3 + 2] = 255;
                    }

                    if (mindistEdge <= 2) {
                        image[((H - y - 1) * W + x) * 3] = 0;
                        image[((H - y - 1) * W + x) * 3 + 1] = 0;
                        image[((H - y - 1) * W + x) * 3 + 2] = 0;
                    }
                }
            }
        }
    }
    std::ostringstream os;
    os << filename << std::setfill('0') << std::setw(4) << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}


void simulate_fluid() {
    int N_fluid = SimulationParameters::N_FLUID;
    int N_air = SimulationParameters::N_AIR;
    int total_particles = N_fluid + N_air;
    std::vector<Particle> particles(total_particles);
    OptimalTransport ot;
    
    ot.diagram.points.resize(total_particles);
    ot.diagram.weights.resize(total_particles);
    ot.diagram.types.resize(total_particles);
    ot.diagram.num_fluid = N_fluid;
    ot.num_fluid_particles = N_fluid;
    ot.num_air_particles = N_air;

    Vector fluid_center(0.5, 0.5, 0.0);
    double desired_fluid_volume = M_PI * SimulationParameters::FLUID_RADIUS * SimulationParameters::FLUID_RADIUS; // different from 0.6 of the course 
    double mass_per_particle = desired_fluid_volume / N_fluid;
    ot.target_fluid_masses.resize(N_fluid, mass_per_particle);
    ot.desired_vol_air = 1.0 - desired_fluid_volume;

    // init fluid particles in circle
    for (int i = 0; i < N_fluid; ++i) {
        particles[i].type = FLUID;
        particles[i].mass = SimulationParameters::MASS;
        particles[i].velocity = Vector(0, 0, 0);
        particles[i].weight = 0.0;

        Vector pos;
        do {
            pos = Vector(uniform(engine), uniform(engine), 0.0);
            pos = (pos - Vector(0.5, 0.5, 0.0)) * 2.0 * SimulationParameters::FLUID_RADIUS + fluid_center;
        } while ((pos - fluid_center).norm2() > SimulationParameters::FLUID_RADIUS * SimulationParameters::FLUID_RADIUS);
        
        particles[i].position = pos;
        particles[i].last_centroid = pos;

        ot.diagram.points[i] = pos;
        ot.diagram.weights[i] = 0.0;
        ot.diagram.types[i] = FLUID;
    }

    // init air particles randomly
    for (int i = 0; i < N_air; ++i) {
        int j = N_fluid + i;
        particles[j].type = AIR;
        particles[j].mass = 0.0;
        particles[j].velocity = Vector(0, 0, 0);
        particles[j].position = Vector(uniform(engine), uniform(engine), 0.0);
        particles[j].weight = 0.0;
        particles[j].last_centroid = particles[j].position;

        ot.diagram.points[j] = particles[j].position;
        ot.diagram.weights[j] = 0.0;
        ot.diagram.types[j] = AIR;
    }
    
    ot.diagram.weights[N_fluid] = 0.0;

    // main loop
    for (int frame = 0; frame < SimulationParameters::NUM_FRAMES; ++frame) {
        std::cout << "Frame " << frame << std::endl;
        for (int i = 0; i < total_particles; ++i) {
            ot.diagram.points[i] = particles[i].position;
            ot.diagram.weights[i] = particles[i].weight;
        }
        ot.optimize();
        for (int i = 0; i < total_particles; ++i) {
            particles[i].weight = ot.diagram.weights[i];
        }
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < N_fluid; ++i) {
            Polygon cell = ot.diagram.cells[i];
            if (cell.vertices.empty()) {
                particles[i].last_centroid = particles[i].position; 
            } else {
                particles[i].last_centroid = cell.centroid();
            }
            
            // apply forces
            Vector spring_force = (1.0 / SimulationParameters::EPSILON_SQ) * (particles[i].last_centroid - particles[i].position);
            Vector gravity_force(0, SimulationParameters::GRAVITY_ACCEL * particles[i].mass, 0);
            Vector total_force = spring_force + gravity_force;
            particles[i].velocity = particles[i].velocity + (SimulationParameters::DT / particles[i].mass) * total_force;
            particles[i].position = particles[i].position + SimulationParameters::DT * particles[i].velocity;

            // we apply here change of direction whenever the particules hits the box
            Vector& pos = particles[i].position;
            Vector& vel = particles[i].velocity;
            if (pos[0] < 0.0) { pos[0] = 0.0; vel[0] *= -1.0; }
            if (pos[0] > 1.0) { pos[0] = 1.0; vel[0] *= -1.0; }
            if (pos[1] < 0.0) { pos[1] = 0.0; vel[1] *= -1.0; }
            if (pos[1] > 1.0) { pos[1] = 1.0; vel[1] *= -1.0; }
        }
        
        save_frame(ot.diagram.cells, particles, "fluid_sim_", frame);
    }
}