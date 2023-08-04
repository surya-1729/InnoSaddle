#ifndef OBJ_PARSER_HPP
#define OBJ_PARSER_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <opencv2/opencv.hpp>
#include <ctime>
#include <filesystem>
#include <string>
/*#include <CGAL/IO/read_off_points.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>*/

namespace fs = std::filesystem;

struct Vec3f {
    float x, y, z;
    Vec3f(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
};

struct Face {
    std::vector<int> vertexIndices;
};

/*
typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef Mesh::Vertex_index Vertex_index;

void printFaces(const Mesh& mesh) {
    for(const auto& face : mesh.faces()) {
        std::cout << "Face: ";
        for(const auto& vertex : vertices_around_face(mesh.halfedge(face), mesh)) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    }
}
*/

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


/*void createMesh(const std::vector<Vec3f>& vertices, const std::vector<Face>& faces, Mesh& mesh) {
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


void readOffFile(const std::string& filePath, std::vector<Vec3f>& vertices, std::vector<Face>& faces) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        return;
    }

    std::string line;
    std::getline(file, line);  // Skip 'OFF' line

    // Read number of vertices, faces, and edges
    std::getline(file, line);
    std::istringstream iss(line);
    int numVertices, numFaces, numEdges;
    iss >> numVertices >> numFaces >> numEdges;

    // Read vertices
    for (int i = 0; i < numVertices; i++) {
        std::getline(file, line);
        std::istringstream iss(line);
        float x, y, z;
        iss >> x >> y >> z;
        vertices.emplace_back(x, y, z);
    }

    // Read faces
    for (int i = 0; i < numFaces; i++) {
        std::getline(file, line);
        std::istringstream iss(line);
        int numVerticesInFace;
        iss >> numVerticesInFace;
        Face face;
        for (int j = 0; j < numVerticesInFace; j++) {
            int vertexIndex;
            iss >> vertexIndex;
            face.vertexIndices.push_back(vertexIndex);
        }
        faces.push_back(face);
    }
}

*/


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
        
    float projectedX = ((vertex.x - left) / (right - left)*width);
    float projectedZ = ((vertex.z - bottom) / (top - bottom)*height);
    float depth = vertex.y; // Use the original z-coordinate as the depth value
    return Vec3f(projectedX, depth, projectedZ);
   
}

bool in_trianglef(double PX, double PZ, double QX, double QZ,
                  double RX, double RZ, double SX, double SZ,
                  float* c1, float* c2)
{
    double z1, z2, n1, n2;

    z1 = (RX * (SZ - PZ) + PX * (RZ - SZ) + SX * (PZ - RZ));
    n1 = (RX * (QZ - PZ) + PX * (RZ - QZ) + QX * (PZ - RZ));

    if ((z1 > 0.0 && n1 < 0.0) || (z1 < 0.0 && n1 > 0.0))
        return false;

    z2 = (QX * (SZ - PZ) + PX * (QZ - SZ) + SX * (PZ - QZ));
    n2 = (QX * (RZ - PZ) + PX * (QZ - RZ) + RX * (PZ - QZ));

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

/*// Function to find the non-zero minimum value in the depth values array
float findNonZeroMin(const std::vector<std::vector<float>>& depthValues) {
    float minDepth = std::numeric_limits<float>::max();

    for (const auto& row : depthValues) {
        for (float value : row) {
            if (value != 0 && value < minDepth) {
                minDepth = value;
            }
        }
    }

    return minDepth;
}*/

// Function to normalize the depth values array
std::vector<std::vector<float>> normalizeDepthValues(const std::vector<std::vector<float>>& depthValues) {
    int width = depthValues.size();
    int height = depthValues[0].size();

    float minDepth = std::numeric_limits<float>::max();
    float maxDepth = std::numeric_limits<float>::lowest();

    for (const auto& row : depthValues) {
        for (float value : row) {
            if (value != 0 && value < minDepth) {
                minDepth = value;
            }
            if (value != 0 && value > maxDepth) {
                maxDepth = value;
            }
        }
    }
    
    //std::cout << "Minimum depth value: " << minDepth << std::endl;
    //std::cout << "Maximum depth value: " << maxDepth << std::endl;

    std::vector<std::vector<float>> normalizedDepthValues(width, std::vector<float>(height));

    for (int x = 0; x < width; ++x) {
        for (int z = 0; z < height; ++z) {
            if (depthValues[x][z] != 0) {
                normalizedDepthValues[x][z] = (depthValues[x][z] - minDepth) / (maxDepth - minDepth);
            } else {
                normalizedDepthValues[x][z] = 0;
            }
        }
    }

    return normalizedDepthValues;
}


std::vector<std::pair<int, uchar>> getMaxPixelValuesAndIndicesInRows(const std::string& imagePath) {
    // Load the image
    cv::Mat img = cv::imread(imagePath, cv::IMREAD_GRAYSCALE);

    std::vector<std::pair<int, uchar>> maxPixelValuesAndIndices;

    for(int i = 0; i < img.rows; i++) {
        double minVal, maxVal;
        cv::Point minLoc, maxLoc;
        cv::minMaxLoc(img.row(i), &minVal, &maxVal, &minLoc, &maxLoc);
        maxPixelValuesAndIndices.push_back({maxLoc.x, static_cast<uchar>(maxVal)});
    }
    return maxPixelValuesAndIndices;
}





void plotPixelImage(const std::vector<std::vector<float>>& normalizedPixelValues, const std::vector<int>& maxDepthIndices, const std::string& fileNameWithoutExtension) {
    



    // Image dimensions
    int width = normalizedPixelValues.size();
    int height = normalizedPixelValues[0].size();

    // Create an OpenCV Mat to store the grayscale pixel image
    cv::Mat grayscaleImage(height, width, CV_8UC1);

    // Generate the pixel image
    for (int x = 0; x < width; ++x) {
        for (int z = 0; z < height; ++z) {
            grayscaleImage.at<uchar>(z, x) = static_cast<uchar>(normalizedPixelValues[x][z]);
        }
    }

    // Create a color image by converting the grayscale image
    cv::Mat pixelImage;
    cv::cvtColor(grayscaleImage, pixelImage, cv::COLOR_GRAY2BGR);

    // Draw the red line with some thickness by connecting adjacent points
    for (int z = 0; z < width - 1; ++z) {
        cv::line(pixelImage, cv::Point(z, maxDepthIndices[z]), cv::Point(z + 1, maxDepthIndices[z + 1]), cv::Scalar(0, 0, 255), 3);
    }

    // Create Mat to store the ridge line
    cv::Mat ridgeLineImage(height, width, CV_8UC1, cv::Scalar(0));
    for (int z = 0; z < width; ++z) {
        ridgeLineImage.at<uchar>(maxDepthIndices[z], z) = 255;
    }

    // Apply Gaussian blur to the ridge line image
    cv::GaussianBlur(ridgeLineImage, ridgeLineImage, cv::Size(5, 5), 2.0, 2.0);

    // Threshold the blurred image to make it binary again
    cv::threshold(ridgeLineImage, ridgeLineImage, 127, 255, cv::THRESH_BINARY);


    // Draw the ridge line with some thickness by connecting adjacent points
    for (int z = 0; z < width - 1; ++z) {
        cv::line(ridgeLineImage, cv::Point(z, maxDepthIndices[z]), cv::Point(z + 1, maxDepthIndices[z + 1]), cv::Scalar(255), 3);
    }

   
    // Generate the file name with date and time
    std::time_t t = std::time(nullptr);
    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "_%Y-%m-%d_%H-%M-%S", std::localtime(&t));
    std::string dateTime = buffer;

    // Create a folder for saving images

    // Create the folder name based on the file name and date-time
    std::string folderName = fileNameWithoutExtension + dateTime;
    std::string folderPath = "pictures/" + folderName;
    fs::create_directory(folderPath);

    // Save the original pixel image
    cv::imwrite(folderPath + "/" + fileNameWithoutExtension + ".png", pixelImage);




    for(int z = 0; z < width; ++z) std::cout << "Pixel value at z=" << z << " and maxDepthIndices[" << z << "]=" << maxDepthIndices[z] << ": " << normalizedPixelValues[z][maxDepthIndices[z]] << std::endl;
    
    std::string filePath = folderPath + "/" + fileNameWithoutExtension + "_values.txt";
    std::ofstream outFile(filePath);

    for(int z = 0; z < width; ++z) {
        outFile << "Pixel value at z=" << z << " and maxDepthIndices[" << z << "]=" << maxDepthIndices[z] << ": " << normalizedPixelValues[z][maxDepthIndices[z]] << std::endl;
    }

    outFile.close();


         // Create a blank image for the profile visualization
    cv::Mat profileImage(height, width, CV_8UC3, cv::Scalar(255, 255, 255)); // White background

    // Determine the maximum and minimum values from normalizedPixelValues
    double maxPixelValue = -std::numeric_limits<double>::infinity();
    double minPixelValue = std::numeric_limits<double>::infinity();
    for (int z = 0; z < width; ++z) {
        double value = normalizedPixelValues[z][maxDepthIndices[z]];
        if (value > maxPixelValue) maxPixelValue = value;
        if (value < minPixelValue) minPixelValue = value;
    }
    double rangePixelValue = maxPixelValue - minPixelValue;

    double scalingFactor = 0.8; // You can adjust this value as needed

    // Draw the profile line
    for (int z = 0; z < width - 1; ++z) {
        int startY = height - 1 - ((normalizedPixelValues[z][maxDepthIndices[z]] - minPixelValue) * height * scalingFactor) / rangePixelValue;
        int endY = height - 1 - ((normalizedPixelValues[z + 1][maxDepthIndices[z + 1]] - minPixelValue) * height * scalingFactor) / rangePixelValue;
        cv::line(profileImage, cv::Point(z, startY), cv::Point(z + 1, endY), cv::Scalar(0, 0, 0), 1); // Black line
    }




    // Save the profile visualization
    std::string profilePath = folderPath + "/" + fileNameWithoutExtension + "_profile.png";
    cv::imwrite(profilePath, profileImage);



    // Increase padding for the graph
    int graphPaddingTop = 50;
    int graphPaddingBottom = 100;
    int graphPaddingLeft = 50;
    int graphPaddingRight = 50;

    // Create graph visualization with additional space for axis labels
    cv::Mat graphImage(height + graphPaddingTop + graphPaddingBottom, width + graphPaddingLeft + graphPaddingRight, CV_8UC3, cv::Scalar(255, 255, 255)); // White background

    // Draw the profile line
    for (int z = 0; z < width - 1; ++z) {
        int startY = graphPaddingTop + height - 1 - ((normalizedPixelValues[z][maxDepthIndices[z]] - minPixelValue) * height * scalingFactor) / rangePixelValue;
        int endY = graphPaddingTop + height - 1 - ((normalizedPixelValues[z + 1][maxDepthIndices[z + 1]] - minPixelValue) * height * scalingFactor) / rangePixelValue;
        cv::line(graphImage, cv::Point(graphPaddingLeft + z, startY), cv::Point(graphPaddingLeft + z + 1, endY), cv::Scalar(0, 0, 0), 1); // Black line
    }

    // Draw axis lines
    cv::line(graphImage, cv::Point(graphPaddingLeft, graphPaddingTop + height), cv::Point(graphPaddingLeft + width, graphPaddingTop + height), cv::Scalar(0, 0, 0), 2);
    cv::line(graphImage, cv::Point(graphPaddingLeft, graphPaddingTop), cv::Point(graphPaddingLeft, graphPaddingTop + height), cv::Scalar(0, 0, 0), 2);

    // Add tick marks and grid lines (example for 10 ticks)
    int numTicks = 10;
    for (int i = 0; i <= numTicks; ++i) {
        int x = graphPaddingLeft + (width * i) / numTicks;
        int y = graphPaddingTop + (height * i) / numTicks;

        // Vertical grid line
        cv::line(graphImage, cv::Point(x, graphPaddingTop), cv::Point(x, graphPaddingTop + height), cv::Scalar(200, 200, 200), 1);

        // Horizontal grid line
        cv::line(graphImage, cv::Point(graphPaddingLeft, y), cv::Point(graphPaddingLeft + width, y), cv::Scalar(200, 200, 200), 1);
    }

    // Add axis labels
    cv::putText(graphImage, "z (index)", cv::Point(graphPaddingLeft + width / 2, graphPaddingTop + height + 50), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(0, 0, 0), 2);
    cv::putText(graphImage, "y depth", cv::Point(10, graphPaddingTop + height / 2), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(0, 0, 0), 2);

    // Save the graph visualization
    std::string graphPath = folderPath + "/" + fileNameWithoutExtension + "_graph.png";
    cv::imwrite(graphPath, graphImage);

    // ...



}













