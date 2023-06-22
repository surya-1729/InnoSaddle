#ifndef OBJ_PARSER_HPP
#define OBJ_PARSER_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>

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

#endif  // OBJ_PARSER_HPP


// MinMax function to get minimum and maximum values of a vector
std::pair<float, float> findMinMax(const std::vector<float>& n) {
    if (n.empty()) {
        // Handle the case when the vector is empty
        std::cerr << "Empty vector." << std::endl;
        return {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
    }

    float minVal = n[0]; // Initialize minVal with the first element of n
    float maxVal = n[0]; // Initialize maxVal with the first element of n

    for (const auto& val : n) {
        if (val < minVal) {
            minVal = val; // Update minVal if a smaller value is found
        }

        if (val > maxVal) {
            maxVal = val; // Update maxVal if a larger value is found
        }
    }

    return {minVal, maxVal};
}

// Orthographic projection function
Vec3f orthographicProjectionf(const Vec3f& vertex, float left, float right, float bottom, float top, int width, int height) {
        
    float projectedX = ((vertex.x - left) / (right - left)*width) - 1;
    float projectedY = ((vertex.y - bottom) / (top - bottom)*height) - 1;
    float depth = vertex.z; // Use the original z-coordinate as the depth value
    return Vec3f(projectedX, projectedY, depth);
   
}

bool in_trianglef(double PX, double PY, double QX, double QY,
                  double RX, double RY, double SX, double SY,
                  float* c1, float* c2)
{
    double z1, z2, n1, n2;

    z1 = (RX * (SY - PY) + PX * (RY - SY) + SX * (PY - RY));
    n1 = (RX * (QY - PY) + PX * (RY - QY) + QX * (PY - RY));

    if ((z1 > 0.0 && n1 < 0.0) || (z1 < 0.0 && n1 > 0.0))
        return false;

    z2 = (QX * (SY - PY) + PX * (QY - SY) + SX * (PY - QY));
    n2 = (QX * (RY - PY) + PX * (QY - RY) + RX * (PY - QY));

    if ((z2 > 0.0 && n2 < 0.0) || (z2 < 0.0 && n2 > 0.0))
        return false;

    if (n1 == 0.0 || n2 == 0.0)
        return false;

    if ((n1 > 0.0) == (n2 > 0.0)) {
        if (z1 * n2 + z2 * n1 <= n1 * n2) {
            *c1 = static_cast<float>(z1 / n1);
            *c2 = static_cast<float>(z2 / n2);
            return true;
        } else {
            return false;
        }
    } else {
        if (z1 * n2 + z2 * n1 >= n1 * n2) {
            *c1 = static_cast<float>(z1 / n1);
            *c2 = static_cast<float>(z2 / n2);
            return true;
        } else {
            return false;
        }
    }
}

// The following is used for calculating the depth value of the pixel
// from the depth values of the corners, using the barycentric coordinates.
float interpolatef(float depthP, float depthQ, float depthR, float c1, float c2)
{
    // Interpolate the depth value using the barycentric coordinates
    float interpolatedDepth = depthP + c1 * (depthQ - depthP) + c2 * (depthR - depthP);
    return interpolatedDepth;
}

