    // Save Sobel derivatives
    cv::imwrite(folderPath + "/" + fileNameWithoutExtension + "_sobelX.png", sobelX);
    cv::imwrite(folderPath + "/" + fileNameWithoutExtension + "_sobelZ.png", sobelZ);
    cv::imwrite(folderPath + "/" + fileNameWithoutExtension + "_sobelXZ.png", sobelXZ);

    // Compute the second-order derivatives in the X and Z directions
    cv::Mat secondOrderX, secondOrderZ;
    cv::Sobel(sobelX, secondOrderX, CV_64F, 1, 0);
    cv::Sobel(sobelZ, secondOrderZ, CV_64F, 0, 1);

    // Analyze the second-order derivatives to find local maximum and minimum
    int localMaxIndex = -1, localMinIndex = -1;
    double localMaxValue = -DBL_MAX, localMinValue = DBL_MAX;

    for (int z = 1; z < width - 1; ++z) {
        // You can choose to use either X or Z direction, depending on the orientation
        double prevValue = secondOrderX.at<double>(maxDepthIndices[z - 1], z - 1);
        double currentValue = secondOrderX.at<double>(maxDepthIndices[z], z);
        double nextValue = secondOrderX.at<double>(maxDepthIndices[z + 1], z + 1);

        // Check for zero crossing from positive to negative (local maximum)
        if (prevValue > 0 && currentValue < 0) {
            if (localMaxValue < prevValue) {
                localMaxValue = prevValue;
                localMaxIndex = z - 1;
            }
        }

        // Check for zero crossing from negative to positive (local minimum)
        if (prevValue < 0 && currentValue > 0) {
            if (localMinValue > prevValue) {
                localMinValue = prevValue;
                localMinIndex = z - 1;
            }
        }
    }

    cv::circle(pixelImage, cv::Point(localMaxIndex, maxDepthIndices[localMaxIndex]), 10, cv::Scalar(0, 255, 0), -1); // green for maximum
    cv::circle(pixelImage, cv::Point(localMinIndex, maxDepthIndices[localMinIndex]), 10, cv::Scalar(255, 0, 0), -1); // blue for minimum

    // Save the modified pixel image
    cv::imwrite(folderPath + "/" + fileNameWithoutExtension + "_localMaxMin.png", pixelImage);

    for(int z = 0; z < width; ++z) std::cout << "Pixel value at z=" << z << " and maxDepthIndices[" << z << "]=" << maxDepthIndices[z] << ": " << normalizedPixelValues[z][maxDepthIndices[z]] << std::endl;
    
    std::string filePath = folderPath + "/" + fileNameWithoutExtension + "_values.txt";
    std::ofstream outFile(filePath);

    for(int z = 0; z < width; ++z) {
        outFile << "Pixel value at z=" << z << " and maxDepthIndices[" << z << "]=" << maxDepthIndices[z] << ": " << normalizedPixelValues[z][maxDepthIndices[z]] << std::endl;
    }

    outFile.close();

    // Create a visualization image
    cv::Mat visualizationImage(height, width, CV_8UC3, cv::Scalar(0, 0, 0)); // Black background

    // Draw the values
    for (int z = 0; z < width; ++z) {
        uchar pixelValue = static_cast<uchar>(normalizedPixelValues[z][maxDepthIndices[z]]);
        cv::line(visualizationImage, cv::Point(z, maxDepthIndices[z]), cv::Point(z, maxDepthIndices[z]), cv::Scalar(pixelValue, pixelValue, pixelValue), 3); // Drawing as grayscale
    }

    // Save the visualization image
    std::string visualizationPath = folderPath + "/" + fileNameWithoutExtension + "_visualization.png";
    cv::imwrite(visualizationPath, visualizationImage);

    // Create a blank image for the profile visualization
    cv::Mat profileImage(height, width, CV_8UC3, cv::Scalar(255, 255, 255)); // White background

    // Normalize the ridge line for proper visualization
    int maxRidgeDepth = *std::max_element(maxDepthIndices.begin(), maxDepthIndices.end());
    int minRidgeDepth = *std::min_element(maxDepthIndices.begin(), maxDepthIndices.end());
    int rangeRidgeDepth = maxRidgeDepth - minRidgeDepth;

    // Draw the profile line
    for (int z = 0; z < width - 1; ++z) {
        int startY = height - 1 - ((maxDepthIndices[z] - minRidgeDepth) * height) / rangeRidgeDepth;
        int endY = height - 1 - ((maxDepthIndices[z + 1] - minRidgeDepth) * height) / rangeRidgeDepth;
        cv::line(profileImage, cv::Point(z, startY), cv::Point(z + 1, endY), cv::Scalar(0, 0, 0), 1); // Black line
    }

    // Save the profile visualization
    std::string profilePath = folderPath + "/" + fileNameWithoutExtension + "_profile.png";
    cv::imwrite(profilePath, profileImage);

    // Create a blank white image
    cv::Mat profileGraph(height, width, CV_8UC3, cv::Scalar(255, 255, 255));

    // Draw the axes
    int axisOffset = 50; // Offset for the axes to leave space for labels
    cv::line(profileGraph, cv::Point(axisOffset, height - axisOffset), cv::Point(width - axisOffset, height - axisOffset), cv::Scalar(0, 0, 0), 2);
    cv::line(profileGraph, cv::Point(axisOffset, height - axisOffset), cv::Point(axisOffset, axisOffset), cv::Scalar(0, 0, 0), 2);

    // Draw the profile line
    for (int z = 0; z < width - 1 - 2 * axisOffset; ++z) {
        int startY = (height - 2 * axisOffset) - 1 - ((maxDepthIndices[z] - minRidgeDepth) * (height - 2 * axisOffset)) / rangeRidgeDepth + axisOffset;
        int endY = (height - 2 * axisOffset) - 1 - ((maxDepthIndices[z + 1] - minRidgeDepth) * (height - 2 * axisOffset)) / rangeRidgeDepth + axisOffset;
        cv::line(profileGraph, cv::Point(z + axisOffset, startY), cv::Point(z + 1 + axisOffset, endY), cv::Scalar(0, 0, 255), 2); // Blue line
    }

    // Optional: Add labels and grid
    // ...


    // ...

    // Draw the axes labels
    int labelOffset = 15; // Offset for labels
    cv::putText(profileGraph, "Y (Depth)", cv::Point(10, axisOffset - labelOffset), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0));
    cv::putText(profileGraph, "X (Index)", cv::Point(width - axisOffset + labelOffset, height - 15), cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 0, 0));

    // Draw tick marks on X and Y axes
    int tickLength = 5; // Length of tick marks
    for (int i = axisOffset; i < width - axisOffset; i += (width - 2 * axisOffset) / 10) {
        cv::line(profileGraph, cv::Point(i, height - axisOffset - tickLength), cv::Point(i, height - axisOffset + tickLength), cv::Scalar(0, 0, 0), 1);
    }
    for (int i = axisOffset; i < height - axisOffset; i += (height - 2 * axisOffset) / 10) {
        cv::line(profileGraph, cv::Point(axisOffset - tickLength, i), cv::Point(axisOffset + tickLength, i), cv::Scalar(0, 0, 0), 1);
    }

    // Optional: Add grid
    // ...

    // Save the plot image
    std::string profileGraphPath = folderPath + "/" + fileNameWithoutExtension + "_profileGraph.png";
    cv::imwrite(profileGraphPath, profileGraph);


    // Find min and max of maxDepthIndices to normalize the plot
    int minDepthValue = *std::min_element(maxDepthIndices.begin(), maxDepthIndices.end());
    int maxDepthValue = *std::max_element(maxDepthIndices.begin(), maxDepthIndices.end());

    // Create an image to hold the plot
    cv::Mat plotImage(height, width, CV_8UC3, cv::Scalar(255, 255, 255)); // white background

    // Plot the points
    for(int x = 0; x < width - 1; ++x) {
        int y1 = height - 1 - ((maxDepthIndices[x] - minDepthValue) * (height - 1) / (maxDepthValue - minDepthValue));
        int y2 = height - 1 - ((maxDepthIndices[x + 1] - minDepthValue) * (height - 1) / (maxDepthValue - minDepthValue));
        cv::line(plotImage, cv::Point(x, y1), cv::Point(x + 1, y2), cv::Scalar(255, 0, 0), 2); // red line
    }

    // Save the plot image
    std::string plotImagePath = folderPath + "/" + fileNameWithoutExtension + "_profileLine.png";
    cv::imwrite(plotImagePath, plotImage);