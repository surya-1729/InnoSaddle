#include <iostream>
#include <vector>
#include "helper.hpp"
#include <limits>

int main() {

    // Path to the OBJ file
    std::string filePath = "./data/processed_data/0002_letti.obj";

    // Variables to store vertices and faces
    std::vector<Vec3f> vertices;
    std::vector<Face> faces;

    // Read the OBJ file
    readObjFile(filePath, vertices, faces);

    //// Parameters for orthographic projection

    // Extract the first column of the matrix
    std::vector<float> firstColumn;
    for (const auto& row : vertices) {
        firstColumn.push_back(row.x);
    }

    // Find the minimum and maximum of the first column
    std::pair<float, float> minMax = findMinMax(firstColumn);
    float left = minMax.first;
    float right = minMax.second;

    // Extract the second column of the matrix
    std::vector<float> secondColumn;
    for (const auto& row : vertices) {
        secondColumn.push_back(row.y);
    }

    // Find the minimum and maximum of the second column
    std::pair<float, float> secondColumnMinMax = findMinMax(secondColumn);
    float bottom = secondColumnMinMax.first;
    float top = secondColumnMinMax.second;

    int width = 1000;
    int height = 1000;

    std::vector<Vec3f> projectedVertices;
    for (const auto& vertex : vertices) {
        Vec3f projectedVertex = orthographicProjectionf(vertex, left, right, bottom, top, width, height);
        projectedVertices.push_back(projectedVertex);
    }
    
    // Initialize a 2D array to store the depth values

    std::vector<std::vector<float>> depthValues(width, std::vector<float>(height, 0.0));
    


    // Iterate over each face
    for (const auto& face : faces) {
        // Calculate the boundaries of the face (min x, max x, min y, max y)
        float minX = std::numeric_limits<float>::max();
        float maxX = std::numeric_limits<float>::lowest();
        float minY = std::numeric_limits<float>::max();
        float maxY = std::numeric_limits<float>::lowest();

        for (const auto& vertexIndex : face.vertexIndices) {
            // Access the projected vertex using the vertex index
            const Vec3f& projectedVertex = projectedVertices[vertexIndex];

            // Update the boundaries based on the projected vertex coordinates
            minX = std::min(minX, projectedVertex.x);
            maxX = std::max(maxX, projectedVertex.x);
            minY = std::min(minY, projectedVertex.y);
            maxY = std::max(maxY, projectedVertex.y);
        }

        // std::cout << minX << " " << maxX << " "  << minY << " " << maxY << "  " << std::endl;

        int minXIndex = static_cast<int>(minX);
        int maxXIndex = static_cast<int>(maxX); 

        int minYIndex = static_cast<int>(minY);
        int maxYIndex = static_cast<int>(maxY); 

        for (int x = minXIndex; x < maxXIndex; ++x) {
            for (int y = minYIndex; y < maxYIndex; ++y) {
                // Check if the pixel center lies inside the triangle
                float c1, c2;
                if (in_trianglef(projectedVertices[face.vertexIndices[0]].x, projectedVertices[face.vertexIndices[0]].y,
                                projectedVertices[face.vertexIndices[1]].x, projectedVertices[face.vertexIndices[1]].y,
                                projectedVertices[face.vertexIndices[2]].x, projectedVertices[face.vertexIndices[2]].y,
                                x + 0.5f, y + 0.5f, &c1, &c2)) {
                    // Pixel center is inside the triangle
                    // Interpolate the depth value using the barycentric coordinates
                    //std::cout << c1 << " 2 " << std::endl;
                    float interpolatedDepth = interpolatef(projectedVertices[face.vertexIndices[0]].z,
                                                        projectedVertices[face.vertexIndices[1]].z,
                                                        projectedVertices[face.vertexIndices[2]].z,
                                                        c1, c2);
                    // Store the depth value in the depth values array
                    
                    depthValues[x][y] = interpolatedDepth;
                    std::cout << "Depth at (" << x << ", " << y << "): " << depthValues[x][y] << std::endl;
                }
            }
                    
        }
 
    }


    return 0;
}
