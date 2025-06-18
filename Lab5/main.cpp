#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cstdio>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "../stb_image.h"

struct Proj {
    double val;
    int idx;
    bool operator<(const Proj& other) const { return val < other.val; }
};

void random_direction(double& x, double& y, double& z) {
    double r1 = (double)rand() / RAND_MAX;
    double r2 = (double)rand() / RAND_MAX;
    x = cos(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
    y = sin(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
    z = 1 - 2 * r2;
}

void sliced_transport(std::vector<double>& input, std::vector<double>& model, int w, int h) {
    int n = w * h;
    
    for (int iter = 0; iter < 50; iter++) {
        double x, y, z;
        random_direction(x, y, z);
        
        std::vector<Proj> proj_in(n), proj_mod(n);
        
        for (int i = 0; i < n; i++) {
            proj_in[i].val = input[i*3] * x + input[i*3+1] * y + input[i*3+2] * z;
            proj_in[i].idx = i;
            proj_mod[i].val = model[i*3] * x + model[i*3+1] * y + model[i*3+2] * z;
            proj_mod[i].idx = i;
        }
        
        std::sort(proj_in.begin(), proj_in.end());
        std::sort(proj_mod.begin(), proj_mod.end());
        
        for (int i = 0; i < n; i++) {
            int idx = proj_in[i].idx;
            double diff = proj_mod[i].val - proj_in[i].val;
            input[idx*3] += diff * x;
            input[idx*3+1] += diff * y;
            input[idx*3+2] += diff * z;
        }
    }
}

int main() {
    int w, h, c, mw, mh, mc;
    unsigned char* img = stbi_load("input.png", &w, &h, &c, STBI_rgb);
    unsigned char* model_img = stbi_load("model.png", &mw, &mh, &mc, STBI_rgb);
    int size = w * h;
    std::vector<double> input(size * 3), model(size * 3);
    for (int i = 0; i < size * 3; i++) {
        input[i] = img[i];
        model[i] = model_img[i];
    }
    
    sliced_transport(input, model, w, h);
    
    std::vector<unsigned char> result(size * 3);
    for (int i = 0; i < size * 3; i++) {
        result[i] = std::max(0.0, std::min(255.0, input[i]));
    }
    
    stbi_write_png("result.png", w, h, 3, &result[0], 0);
    
    return 0;
}