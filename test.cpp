#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <fstream>
#include <vector>
#include <filesystem>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef Mesh::Vertex_index Vertex_index;

struct Vec3f {
    float x, y, z;
    Vec3f(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
};

struct Face {
    std::vector<int> vertexIndices;
};

void readObjFile(const std::string& filePath, std::vector<Vec3f>& vertices, std::vector<Face>& faces) {
    // Open the OBJ file
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        return;
    }

    // Read the file line by line
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if (prefix == "v") {
            // Read vertex coordinates
            float x, y, z;
            iss >> x >> y >> z;
            vertices.emplace_back(x, y, z);
        } else if (prefix == "f") {
            // Read face vertex indices
            std::string vertexStr;
            Face face;
            while (iss >> vertexStr) {
                std::istringstream viss(vertexStr);
                int vertexIndex;
                if (viss >> vertexIndex) {
                    face.vertexIndices.push_back(vertexIndex - 1);  // OBJ indices are 1-based
                }
            }
            faces.push_back(face);
        }
    }
}

void createMesh(const std::vector<Vec3f>& vertices, const std::vector<Face>& faces, Mesh& mesh) {
    std::vector<Vertex_index> vertex_indices;

    // Add vertices to the mesh
    for (const Vec3f& vertex : vertices) {
        vertex_indices.push_back(mesh.add_vertex(K::Point_3(vertex.x, vertex.y, vertex.z)));
    }

    // Add faces to the mesh
    for (const Face& face : faces) {
        std::vector<Vertex_index> face_indices;
        for (int index : face.vertexIndices) {
            face_indices.push_back(vertex_indices[index]);
        }
        mesh.add_face(face_indices);
    }
}

void writeMeshToFile(const Mesh& mesh, const std::string& filePath) {
    std::ofstream file(filePath, std::ios::out);
    if (!file) {
        std::cerr << "Could not open file for writing: " << filePath << std::endl;
        return;
    }

    if (!CGAL::IO::write_OFF(file, mesh)) {
        std::cerr << "Error writing to file: " << filePath << std::endl;
    }
}


int main() {
    std::string filePath = "./data/processed_data/0004_piri.obj";
    std::vector<Vec3f> vertices;
    std::vector<Face> faces;

    // Extract base file name from input file path
    std::filesystem::path inputPath(filePath);
    std::string baseName = inputPath.stem().string();

    readObjFile(filePath, vertices, faces);

    Mesh mesh;
    createMesh(vertices, faces, mesh);

    CGAL::Polygon_mesh_processing::triangulate_faces(mesh);

    // Now the mesh is ready for further processing
    // ...
    // Fill the holes in the mesh.
    unsigned int nb_holes = 0;
    for(auto h : halfedges(mesh)) {
        if(is_border(h, mesh)) {
            CGAL::Euler::fill_hole(h, mesh);
            ++nb_holes;
        }
    }
    std::cout << "Filled " << nb_holes << " holes.\n";
    // Append new extension to base file name to create output file path
    std::string outputFilePath = "/home/surya/cpp_project/Innosaddle/data/updated/" + baseName + "_.off";
    writeMeshToFile(mesh, outputFilePath);
    return 0;
}




