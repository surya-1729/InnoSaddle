#include <iostream>
#include <vector>
#include "helper.hpp"
#include <limits>
#include <filesystem>
#include <string>
// #include <opencv2/opencv.hpp>

int main() {

    // Path to the OBJ file
    std::string filePath = "./data/processed_data/0003_maicon.obj";

    // Extract the file name without extension
    std::filesystem::path path(filePath);
    std::string fileNameWithoutExtension = path.stem().string();

    // Variables to store vertices and faces
    std::vector<Vec3f> vertices;
    std::vector<Face> faces;

    
    // Read the OBJ file
    readObjFile(filePath, vertices, faces);

    /*
    std::cout << "Number of faces: " << createfaces.size() << std::endl;

    Mesh mesh;
    createMesh(createvertices, createfaces, mesh);

    CGAL::Polygon_mesh_processing::triangulate_faces(mesh);
    //CGAL::Polygon_mesh_processing::triangulate_hole_polyline(mesh)

    // Fill the holes in the mesh.
    unsigned int nb_holes = 0;
    for(auto h : halfedges(mesh)) {
        if(is_border(h, mesh)) {
            CGAL::Euler::fill_hole(h, mesh);
            ++nb_holes;
        }
    }
    std::cout << "Filled " << nb_holes << " holes.\n";

    std::cout << "Number of faces after filling holes: " << mesh.num_faces() << std::endl;

    // Triangulate again after filling holes.
    // CGAL::Polygon_mesh_processing::triangulate_faces(mesh)

    unsigned int nb_NOTTRIANGULATED = 0;
    // Iterate over each face in the mesh
    for(auto f : faces(mesh))
    {
        // Triangulate the face
        bool success = CGAL::Polygon_mesh_processing::triangulate_face(f, mesh);

        // Check if triangulation was successful
        if (!success)
        {
            ++nb_NOTTRIANGULATED;
            //std::cerr << "Failed to triangulate face " << f << std::endl;
            //std::cout << "Filled " << nb_NOTTRIANGULATED << " holes.\n";
        }
    }


    printFaces(mesh); // Print the faces

    // Append new extension to base file name to create output file path
    std::string outputFilePath = "/home/surya/cpp_project/Innosaddle/data/processed_data/" + fileNameWithoutExtension + "_.off";
    writeMeshToFile(mesh, outputFilePath);

    // Variables to store vertices and faces
    std::vector<Vec3f> vertices;
    std::vector<Face> faces;

    // Read the OFF file
    readOffFile(outputFilePath, vertices, faces);*/

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

    // Extract the third column of the matrix
    std::vector<float> thirdColumn;
    for (const auto& row : vertices) {
        thirdColumn.push_back(row.z);
    }

    // Find the minimum and maximum of the third column
    std::pair<float, float> thirdColumnMinMax = findMinMax(thirdColumn);
    float bottom = thirdColumnMinMax.first;
    float top = thirdColumnMinMax.second;

    int width = 1000;
    int height = 1000;

    std::vector<Vec3f> projectedVertices;
    for (const auto& vertex : vertices) {
        Vec3f projectedVertex = orthographicProjectionf(vertex, left, right, bottom, top, width, height);
        projectedVertices.push_back(projectedVertex);
    }


    // Initialize a 2D array to store the depth values
    std::vector<std::vector<float>> depthValues(width, std::vector<float>(height, 0.0));
    int printCount = 0; // Add this before you start iterating over the faces
    // Iterate over each face
    for (const auto& face : faces) {
        // Calculate the boundaries of the face (min x, max x, min z, max z)
        float minX = std::numeric_limits<float>::max();
        float maxX = std::numeric_limits<float>::lowest();
        float minZ = std::numeric_limits<float>::max();
        float maxZ = std::numeric_limits<float>::lowest();

        for (const auto& vertexIndex : face.vertexIndices) {
            // Access the projected vertex using the vertex index
            const Vec3f& projectedVertex = projectedVertices[vertexIndex];

            // Update the boundaries based on the projected vertex coordinates
            minX = std::min(minX, projectedVertex.x);
            maxX = std::max(maxX, projectedVertex.x);
            minZ = std::min(minZ, projectedVertex.z);
        maxZ = std::max(maxZ, projectedVertex.z);   
        }


        int minXIndex = static_cast<int>(minX);
        int maxXIndex = static_cast<int>(maxX); 

        int minZIndex = static_cast<int>(minZ);
        int maxZIndex = static_cast<int>(maxZ); 

        for (int x = minXIndex; x < maxXIndex + 1; ++x) {
            for (int z = minZIndex; z < maxZIndex + 1; ++z) {
                // Check if the pixel center lies inside the triangle
                float c1, c2;
                if (in_trianglef(projectedVertices[face.vertexIndices[0]].x, projectedVertices[face.vertexIndices[0]].z,
                                projectedVertices[face.vertexIndices[1]].x, projectedVertices[face.vertexIndices[1]].z,
                                projectedVertices[face.vertexIndices[2]].x, projectedVertices[face.vertexIndices[2]].z,
                                x + 0.5f, z + 0.5f, &c1, &c2)) {
                    // Pixel center is inside the triangle
                    // Interpolate the depth value using the barycentric coordinates
                    //std::cout << c1 << " 2 " << std::endl;
                    float interpolatedDepth = interpolatef(projectedVertices[face.vertexIndices[0]].y,
                                                        projectedVertices[face.vertexIndices[1]].y,
                                                        projectedVertices[face.vertexIndices[2]].y,
                                                        c1, c2);

                    /*// Add the following lines
                    if (printCount < 100000000) { // Only print the first 10
                        std::cout << "Interpolated depth for face " << printCount << ": " << interpolatedDepth << "\n";
                        ++printCount;
                    }*/
                    // Store the depth value in the depth values array                    
                    depthValues[x][z] = interpolatedDepth;
                }
            }
        }
    }



    // Multiply each value in the depth values array by -1
    /*for (auto& row : depthValues) {
        for (float& value : row) {
            value *= -1;
        }
    }*/


    std::vector<std::vector<float>> normalizedDepthValues = normalizeDepthValues(depthValues);
    std::vector<std::vector<float>> normalizedPixelValues = normalizedDepthValues;  // Create a new array for storing the multiplied values

    // Multiply each value by 255 and store in normalizedPixelValues
    for (auto& row : normalizedPixelValues) {
        for (float& value : row) {
            value *= 255;
        }
    }

    /*for (auto& row : normalizedPixelValues) {
    for (float& value : row) {
        value = (value > 100) ? 255 : 0;
    }
    }*/


    /*// After calculating normalizedPixelValues
    for (int i = 0; i < normalizedPixelValues.size(); ++i) {
        for (int j = 0; j < normalizedPixelValues[i].size(); ++j) {
            std::cout << "Normalized pixel value at (" << i << ", " << j << "): " << normalizedPixelValues[i][j] << "\n";
        }
    }*/


    // Initialize an array to store the pixel indices of maximum depth values in each row
    std::vector<int> maxDepthIndices(width, 0);

    // Iterate over each row of the depthValues array
    for (int z = 0; z < width; ++z) {
        float maxDepth = std::numeric_limits<float>::lowest();
        int maxDepthIndex = 0;

        for (int x = 0; x < height; ++x) {
            // Check if the current depth value is greater than the maximum depth
            if (normalizedPixelValues[z][x] > maxDepth) {
                maxDepth = normalizedPixelValues[z][x];
                maxDepthIndex = x;
            }
        }

        // Store the pixel index with maximum depth value for the current row
        maxDepthIndices[z] = maxDepthIndex;
    }

    for (int i = 0; i < maxDepthIndices.size(); ++i) {
    std::cout << "Row " << i << ": Maximum depth at index = " << maxDepthIndices[i] << std::endl;
    }


    plotPixelImage(normalizedPixelValues, maxDepthIndices, fileNameWithoutExtension);
    return 0;
}
