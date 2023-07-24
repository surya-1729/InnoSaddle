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
#include <CGAL/IO/read_off_points.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

namespace fs = std::filesystem;

struct Vec3f {
    float x, y, z;
    Vec3f(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
};

struct Face {
    std::vector<int> vertexIndices;
};


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




void plotPixelImage(const std::vector<std::vector<float>>& normalizedPixelValues, const std::string& fileNameWithoutExtension) {
    // Image dimensions
    int width = normalizedPixelValues.size();
    int height = normalizedPixelValues[0].size();

    // Create an OpenCV Mat to store the pixel image
    cv::Mat pixelImage(height, width, CV_8UC1);

    // Generate the pixel image
    for (int x = 0; x < width; ++x) {
        for (int z = 0; z < height; ++z) {
            // Set the pixel value from the normalized pixel values array
            pixelImage.at<uchar>(z, x) = static_cast<uchar>(normalizedPixelValues[x][z]);
        }
    }

    // Generate the file name with date and time
    std::time_t t = std::time(nullptr);
    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "_%Y-%m-%d_%H-%M-%S", std::localtime(&t));
    std::string dateTime = buffer;

    // Create the folder name based on the file name and date-time
    std::string folderName = fileNameWithoutExtension + dateTime;

    // Create the directory
    std::string folderPath = "pictures/" + folderName;
    fs::create_directory(folderPath);

    // Generate the complete file name for the pixel image
    std::string pixelFileName = fileNameWithoutExtension + ".png";

    // Save the pixel image
    std::string pixelFilePath = folderPath + "/" + pixelFileName;
    cv::imwrite(pixelFilePath, pixelImage);

    // Perform Gaussian blur
    cv::Mat blurredImage;
    cv::GaussianBlur(pixelImage, blurredImage, cv::Size(5, 5), 0);

    // Generate the complete file name for the blurred image
    std::string blurredFileName = fileNameWithoutExtension + "_blurred.png";

    // Save the blurred image
    std::string blurredFilePath = folderPath + "/" + blurredFileName;
    cv::imwrite(blurredFilePath, blurredImage);

    // Perform Contrast increase to blurred Image
    cv::Mat contrastImage;
    cv::equalizeHist(blurredImage, contrastImage);

    // Generate the complete file name for the contrast image
    std::string contrastFileName = fileNameWithoutExtension + "_contrast.png";

    // Save the contrast image
    std::string contrastFilePath = folderPath + "/" + contrastFileName;
    cv::imwrite(contrastFilePath, contrastImage);


    // Perform Canny edge detection
    cv::Mat cannyImage;
    cv::Canny(contrastImage, cannyImage, 50, 150);

    // Generate the complete file name for the canny image
    std::string cannyFileName = fileNameWithoutExtension + "_canny.png";

    // Save the canny image
    std::string cannyFilePath = folderPath + "/" + cannyFileName;
    cv::imwrite(cannyFilePath, cannyImage);

    // Perform Sobel operator in the x-direction
    cv::Mat sobelXImage;
    cv::Sobel(contrastImage, sobelXImage, CV_8U, 1, 0);

    // Generate the complete file name for the sobel x-image
    std::string sobelXFileName = fileNameWithoutExtension + "_sobel_x.png";

    // Save the sobel x-image
    std::string sobelXFilePath = folderPath + "/" + sobelXFileName;
    cv::imwrite(sobelXFilePath, sobelXImage);

    // Perform Sobel operator in the y-direction
    cv::Mat sobelYImage;
    cv::Sobel(contrastImage, sobelYImage, CV_8U, 0, 1);

    // Generate the complete file name for the sobel y-image
    std::string sobelYFileName = fileNameWithoutExtension + "_sobel_y.png";

    // Save the sobel y-image
    std::string sobelYFilePath = folderPath + "/" + sobelYFileName;
    cv::imwrite(sobelYFilePath, sobelYImage);

    // Perform Sobel operator in the xy-direction
    cv::Mat sobelXYImage;
    cv::Sobel(contrastImage, sobelXYImage, CV_8U, 1, 1);

    // Generate the complete file name for the sobel y-image
    std::string sobelXYFileName = fileNameWithoutExtension + "_sobel_xy.png";

    // Save the sobel xy-image
    std::string sobelXYFilePath = folderPath + "/" + sobelXYFileName;
    cv::imwrite(sobelXYFilePath, sobelXYImage);

    // Assign colors based on intensity values for Sobel x-image
    cv::Mat sobelXColorImage;
    cv::applyColorMap(sobelXImage, sobelXColorImage, cv::COLORMAP_HOT);

    // Generate the complete file name for the colored Sobel x-image
    std::string sobelXColorFileName = fileNameWithoutExtension + "_sobel_x_color.png";

    // Save the colored Sobel x-image
    std::string sobelXColorFilePath = folderPath + "/" + sobelXColorFileName;
    cv::imwrite(sobelXColorFilePath, sobelXColorImage);

    // Assign colors based on intensity values for Sobel y-image
    cv::Mat sobelYColorImage;
    cv::applyColorMap(sobelYImage, sobelYColorImage, cv::COLORMAP_HOT);

    // Generate the complete file name for the colored Sobel y-image
    std::string sobelYColorFileName = fileNameWithoutExtension + "_sobel_y_color.png";

    // Save the colored Sobel y-image
    std::string sobelYColorFilePath = folderPath + "/" + sobelYColorFileName;
    cv::imwrite(sobelYColorFilePath, sobelYColorImage);

    // Assign colors based on intensity values for Sobel xy-image
    cv::Mat sobelXYColorImage;
    cv::applyColorMap(sobelXYImage, sobelXYColorImage, cv::COLORMAP_HOT);

    // Generate the complete file name for the colored Sobel xy-image
    std::string sobelXYColorFileName = fileNameWithoutExtension + "_sobel_xy_color.png";

    // Save the colored Sobel xy-image
    std::string sobelXYColorFilePath = folderPath + "/" + sobelXYColorFileName;
    cv::imwrite(sobelXYColorFilePath, sobelXYColorImage);

    // Perform Laplacian operator
    cv::Mat laplacianImage;
    cv::Laplacian(contrastImage, laplacianImage, CV_8U);

    // Generate the complete file name for the laplacian image
    std::string laplacianFileName = fileNameWithoutExtension + "_laplacian.png";

    // Save the laplacian image
    std::string laplacianFilePath = folderPath + "/" + laplacianFileName;
    cv::imwrite(laplacianFilePath, laplacianImage);

    // Determine the maximum width and height among the images
    int maxWidth = std::max(sobelXImage.cols, std::max(sobelYImage.cols, std::max(sobelXYImage.cols, std::max(laplacianImage.cols, std::max(cannyImage.cols, contrastImage.cols)))));
    int maxHeight = std::max(sobelXImage.rows, std::max(sobelYImage.rows, std::max(sobelXYImage.rows, std::max(laplacianImage.rows, std::max(cannyImage.rows, contrastImage.rows)))));


    // Create a canvas to hold the combined images and legends
    int canvasWidth = 3 * maxWidth;
    int canvasHeight = 2 * maxHeight + 100; // Add space for legends
    cv::Mat canvas(canvasHeight, canvasWidth, CV_8UC1, cv::Scalar(255));

    // Copy the Sobel x-image to the canvas
    cv::Mat canvasROI_X = canvas(cv::Rect(0, 0, sobelXImage.cols, sobelXImage.rows));
    sobelXImage.copyTo(canvasROI_X);
    cv::putText(canvas, "Sobel X", cv::Point(450, maxHeight + 30), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(0, 255, 0), 2);

    // Copy the Sobel y-image to the canvas
    cv::Mat canvasROI_Y = canvas(cv::Rect(maxWidth, 0, sobelYImage.cols, sobelYImage.rows));
    sobelYImage.copyTo(canvasROI_Y);
    cv::putText(canvas, "Sobel Y", cv::Point(maxWidth + 450, maxHeight + 30), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(0, 255, 0), 2);

    // Copy the Sobel xy-image to the canvas
    cv::Mat canvasROI_XY = canvas(cv::Rect(2 * maxWidth, 0, sobelXYImage.cols, sobelXYImage.rows));
    sobelXYImage.copyTo(canvasROI_XY);
    cv::putText(canvas, "Sobel XY", cv::Point(2 * maxWidth + 450, maxHeight + 30), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(0, 255, 0), 2);

    // Copy the laplacian image to the canvas
    cv::Mat canvasROI_lap = canvas(cv::Rect(0, maxHeight + 100, laplacianImage.cols, laplacianImage.rows));
    laplacianImage.copyTo(canvasROI_lap);
    cv::putText(canvas, "Laplacian", cv::Point(450, maxHeight + 70), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(0, 255, 0), 2);

    // Copy the Canny image to the canvas
    cv::Mat canvasROI_canny = canvas(cv::Rect(maxWidth, maxHeight + 100, cannyImage.cols, cannyImage.rows));
    cannyImage.copyTo(canvasROI_canny);
    cv::putText(canvas, "Canny", cv::Point(maxWidth + 450, maxHeight + 70), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(0, 255, 0), 2);

    // Copy the Contrast image to the canvas
    cv::Mat canvasROI_contrast = canvas(cv::Rect(2 * maxWidth, maxHeight + 100, contrastImage.cols, contrastImage.rows));
    contrastImage.copyTo(canvasROI_contrast);
    cv::putText(canvas, "Contrast", cv::Point(2 * maxWidth + 450, maxHeight + 70), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(0, 255, 0), 2);

    // Generate the complete file name for the combined image
    std::string combinedFileName = fileNameWithoutExtension + "_combined.png";

    // Save the combined image
    std::string combinedFilePath = folderPath + "/" + combinedFileName;
    cv::imwrite(combinedFilePath, canvas);

}











