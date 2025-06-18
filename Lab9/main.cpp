#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <map>
#include <utility>
#include <sstream>
#include <iomanip>
#include <set>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Vector {
    double x, y, z;

    Vector(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0) : x(x_), y(y_), z(z_) {}

    Vector operator+(const Vector& other) const { return Vector(x + other.x, y + other.y, z + other.z); }
    Vector operator-(const Vector& other) const { return Vector(x - other.x, y - other.y, z - other.z); }
    Vector operator*(double scalar) const { return Vector(x * scalar, y * scalar, z * scalar); }
    Vector operator/(double scalar) const { return Vector(x / scalar, y / scalar, z / scalar); }
    double length() const { return std::sqrt(x * x + y * y + z * z); }
};

struct TriangleIndices {
    int vtxi, vtxj, vtxk;
};

class TriangleMesh {
public:
    std::vector<Vector> vertices;
    std::vector<TriangleIndices> indices;
    std::vector<std::vector<int>> adjacencyList;

    void readOBJ(const char* obj) {
        std::ifstream ifs(obj);
        if (!ifs.is_open()) {
            std::cerr << "Error: Could not open file " << obj << std::endl;
            return;
        }

        std::string line;
        while (std::getline(ifs, line)) {
            std::stringstream ss(line);
            std::string token;
            ss >> token;

            if (token == "v") {
                Vector vec;
                if (ss >> vec.x >> vec.y >> vec.z) {
                    vertices.push_back(vec);
                }
            } else if (token == "f") {
                std::string face_str = line.substr(2);
                std::stringstream face_ss(face_str);
                std::vector<int> face_vertices;
                
                std::string vertex_data;
                while (face_ss >> vertex_data) {
                    int vertex_idx = std::stoi(vertex_data.substr(0, vertex_data.find('/')));
                    face_vertices.push_back(vertex_idx - 1);
                }

                for (size_t i = 1; i < face_vertices.size() - 1; ++i) {
                    TriangleIndices tri;
                    tri.vtxi = face_vertices[0];
                    tri.vtxj = face_vertices[i];
                    tri.vtxk = face_vertices[i + 1];
                    indices.push_back(tri);
                }
            }
        }
        ifs.close();
        buildAdjacencyList();
    }

    void buildAdjacencyList() {
        adjacencyList.clear();
        adjacencyList.resize(vertices.size());
        
        for (const auto& tri : indices) {
            adjacencyList[tri.vtxi].push_back(tri.vtxj);
            adjacencyList[tri.vtxi].push_back(tri.vtxk);
            adjacencyList[tri.vtxj].push_back(tri.vtxi);
            adjacencyList[tri.vtxj].push_back(tri.vtxk);
            adjacencyList[tri.vtxk].push_back(tri.vtxi);
            adjacencyList[tri.vtxk].push_back(tri.vtxj);
        }

        for (auto& neighbors : adjacencyList) {
            std::sort(neighbors.begin(), neighbors.end());
            neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
        }
    }

    std::vector<int> getBoundaryVertices() const {
        std::map<std::pair<int, int>, int> edge_counts;
        
        for (const auto& tri : indices) {
            std::pair<int, int> edges[] = {
                {std::min(tri.vtxi, tri.vtxj), std::max(tri.vtxi, tri.vtxj)},
                {std::min(tri.vtxj, tri.vtxk), std::max(tri.vtxj, tri.vtxk)},
                {std::min(tri.vtxk, tri.vtxi), std::max(tri.vtxk, tri.vtxi)}
            };

            for (const auto& edge : edges) {
                edge_counts[edge]++;
            }
        }

        std::set<int> boundary_vertex_set;
        std::map<int, std::vector<int>> boundary_adjacency;
        
        for (const auto& pair : edge_counts) {
            if (pair.second == 1) {
                int v1 = pair.first.first;
                int v2 = pair.first.second;
                boundary_vertex_set.insert(v1);
                boundary_vertex_set.insert(v2);
                boundary_adjacency[v1].push_back(v2);
                boundary_adjacency[v2].push_back(v1);
            }
        }

        if (boundary_vertex_set.empty()) {
            return {};
        }

        std::vector<int> boundary_vertices;
        int start_vertex = *boundary_vertex_set.begin();
        int current = start_vertex;
        int previous = -1;
        
        do {
            boundary_vertices.push_back(current);
            int next = -1;
            for (int neighbor : boundary_adjacency[current]) {
                if (neighbor != previous) {
                    next = neighbor;
                    break;
                }
            }
            previous = current;
            current = next;
        } while (current != start_vertex && current != -1);

        return boundary_vertices;
    }

    const std::vector<int>& getNeighbors(int vertexIndex) const {
        return adjacencyList[vertexIndex];
    }

    std::vector<Vector> tutteEmbedding(int nbiter = 500, double tolerance = 1e-6) {
        if (vertices.empty()) {
            std::cerr << "Error: Mesh has no vertices." << std::endl;
            return {};
        }

        std::vector<Vector> embeddedVertices = vertices;
        std::vector<int> boundaryVerticesIndices = getBoundaryVertices();
        
        if (boundaryVerticesIndices.empty()) {
            std::cerr << "Error: No boundary vertices found." << std::endl;
            return {};
        }

        int n_boundary = boundaryVerticesIndices.size();
        
        // Calculate total boundary length
        double s = 0.0;
        for (int i = 0; i < n_boundary; ++i) {
            int curr = boundaryVerticesIndices[i];
            int next = boundaryVerticesIndices[(i + 1) % n_boundary];
            s += (vertices[next] - vertices[curr]).length();
        }
        
        double cs = 0.0;
        for (int i = 0; i < n_boundary; ++i) {
            int b_idx = boundaryVerticesIndices[i];
            double theta = 2.0 * M_PI * cs / s;
            embeddedVertices[b_idx].x = std::cos(theta);
            embeddedVertices[b_idx].y = std::sin(theta);
            embeddedVertices[b_idx].z = 0.0;
            int next = boundaryVerticesIndices[(i + 1) % n_boundary];
            cs += (vertices[next] - vertices[b_idx]).length();
        }

        std::vector<bool> isBoundary(vertices.size(), false);
        for (int idx : boundaryVerticesIndices) {
            isBoundary[idx] = true;
        }

        // Jacobi iterations 
        for (int iter = 0; iter < nbiter; ++iter) {
            std::vector<Vector> newPositions = embeddedVertices;
            double max_displacement = 0.0;
            for (size_t i = 0; i < vertices.size(); ++i) {
                if (!isBoundary[i]) {
                    const std::vector<int>& neighbors = getNeighbors(i);
                    if (!neighbors.empty()) {
                        Vector sum(0, 0, 0);
                        for (int neighborIdx : neighbors) {
                            sum = sum + embeddedVertices[neighborIdx]; 
                        }
                        Vector newPos = sum / neighbors.size(); 
                        max_displacement = std::max(max_displacement, (newPos - embeddedVertices[i]).length());
                        newPositions[i] = newPos;
                    }
                }
            }
            embeddedVertices = newPositions;

            if (max_displacement < tolerance) {
                break;
            }
        }

        return embeddedVertices;
    }

    void writeOBJ(const char* filename, const std::vector<Vector>& outputVertices) const {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
            return;
        }
        outFile << std::fixed << std::setprecision(8);
        for (const auto& v : outputVertices) {
            outFile << "v " << v.x << " " << v.y << " " << v.z << "\n";
        }
        for (const auto& tri : indices) {
            outFile << "f " << tri.vtxi + 1 << " " << tri.vtxj + 1 << " " << tri.vtxk + 1 << "\n";
        }
        outFile.close();
    }
};

int main() {
    TriangleMesh mesh;
    mesh.readOBJ("Goethe.obj");

    if (mesh.vertices.empty()) {
        std::cerr << "Failed to load mesh." << std::endl;
        return 1;
    }
    std::vector<Vector> tutteResult = mesh.tutteEmbedding(1000, 1e-7);
    if (tutteResult.empty()) {
        std::cerr << "Tutte's embedding failed." << std::endl;
        return 1;
    }
    mesh.writeOBJ("goethe_tutte_embedding.obj", tutteResult);
    return 0;
}