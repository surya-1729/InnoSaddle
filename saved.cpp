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













    // Create graph visualization with additional space for axis labels
    cv::Mat graphImage(height + 100, width, CV_8UC3, cv::Scalar(255, 255, 255)); // White background


    // Draw the same profile line
    for (int z = 0; z < width - 1; ++z) {
        int startY = height - 1 - ((normalizedPixelValues[z][maxDepthIndices[z]] - minPixelValue) * height * scalingFactor) / rangePixelValue;
        int endY = height - 1 - ((normalizedPixelValues[z + 1][maxDepthIndices[z + 1]] - minPixelValue) * height * scalingFactor) / rangePixelValue;
        cv::line(graphImage, cv::Point(z, startY), cv::Point(z + 1, endY), cv::Scalar(0, 0, 0), 1); // Black line
    }


    // Draw axis lines
    cv::line(graphImage, cv::Point(0, height), cv::Point(width, height), cv::Scalar(0, 0, 0), 2);
    cv::line(graphImage, cv::Point(0, 0), cv::Point(0, height), cv::Scalar(0, 0, 0), 2);

    // Add axis labels
    cv::putText(graphImage, "z (index)", cv::Point(width / 2, height + 50), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(0, 0, 0), 2);
    cv::putText(graphImage, "y depth", cv::Point(10, height / 2), cv::FONT_HERSHEY_SIMPLEX, 0.7, cv::Scalar(0, 0, 0), 2);

    // Save the graph visualization
    std::string graphPath = folderPath + "/" + fileNameWithoutExtension + "_graph.png";
    cv::imwrite(graphPath, graphImage);






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






     std::vector<double> firstDerivative(width - 1);
    std::vector<double> secondDerivative(width - 2);

    // Compute first and second derivatives
    for (int z = 1; z < width; ++z) {
        firstDerivative[z - 1] = profile[z] - profile[z - 1];
    }
    for (int z = 1; z < width - 1; ++z) {
        secondDerivative[z - 1] = firstDerivative[z] - firstDerivative[z - 1];
    }

    // Find and draw local maxima and minima
    for (int z = 1; z < width - 1; ++z) {
        if (firstDerivative[z - 1] * firstDerivative[z] < 0) { // Sign change in first derivative
            if (secondDerivative[z - 1] > 0) {
                // Local minimum
                int y = height - 1 - ((profile[z] - minPixelValue) * height * scalingFactor) / rangePixelValue;
                cv::circle(graphImage, cv::Point(z, y), 5, cv::Scalar(255, 0, 0), -1); // Draw red circle
            } else if (secondDerivative[z - 1] < 0) {
                // Local maximum
                int y = height - 1 - ((profile[z] - minPixelValue) * height * scalingFactor) / rangePixelValue;
                cv::circle(graphImage, cv::Point(z, y), 5, cv::Scalar(0, 255, 0), -1); // Draw green circle
            }
        }
    }


    // Save the graph visualization
    std::string graphPath = folderPath + "/" + fileNameWithoutExtension + "_graph.png";
    cv::imwrite(graphPath, graphImage);





    // Apply Gaussian smoothing to the profile
    cv::Mat profileMat(profile.size(), 1, CV_64F, profile.data());
    cv::GaussianBlur(profileMat, profileMat, cv::Size(5, 1), 2.0, 2.0);

    // Compute first and second derivatives
    for (int z = 1; z < width; ++z) {
        firstDerivative[z - 1] = profile[z] - profile[z - 1];
    }
    for (int z = 1; z < width - 1; ++z) {
        secondDerivative[z - 1] = firstDerivative[z] - firstDerivative[z - 1];
    }

    // Find the global maximum and minimum of the second derivative
    int maxIdx = std::distance(secondDerivative.begin(), std::max_element(secondDerivative.begin(), secondDerivative.end()));
    int minIdx = std::distance(secondDerivative.begin(), std::min_element(secondDerivative.begin(), secondDerivative.end()));

    // The global maxima and minima of the second derivative correspond to the local maximum and minimum of the original profile
    int yMax = height - 1 - ((profile[maxIdx + 1] - minPixelValue) * height * scalingFactor) / rangePixelValue;
    int yMin = height - 1 - ((profile[minIdx + 1] - minPixelValue) * height * scalingFactor) / rangePixelValue;
    cv::circle(graphImage, cv::Point(graphPaddingLeft + maxIdx + 1, yMax), 5, cv::Scalar(0, 255, 0), -1); // Draw green circle for local max
    cv::circle(graphImage, cv::Point(graphPaddingLeft + minIdx + 1, yMin), 5, cv::Scalar(255, 0, 0), -1); // Draw red circle for local min

    // Save the graph visualization
    std::string graphPath = folderPath + "/" + fileNameWithoutExtension + "_graph.png";
    cv::imwrite(graphPath, graphImage);






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















 







 std::vector<double> profile(width);
    for (int z = 0; z < width; ++z) {
        profile[z] = normalizedPixelValues[z][maxDepthIndices[z]];
    }

   
    // Degree of polynomial
    int degree = 5;

    // Design matrix
    Eigen::MatrixXd X(profile.size(), degree + 1);
    for (int i = 0; i < profile.size(); ++i) {
        for (int j = 0; j <= degree; ++j) {
            X(i, j) = std::pow(i, j);
        }
    }

    // Profile as a vector
    Eigen::VectorXd y(profile.size());
    for (int i = 0; i < profile.size(); ++i) {
        y(i) = profile[i];
    }

    // Solve for coefficients
    Eigen::VectorXd coefficients = X.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);

    // Compute the first and second derivatives at a given point z
    double z = 10; // Example point
    double y_prime = 0;
    double y_double_prime = 0;
    for (int j = 1; j <= degree; ++j) {
        y_prime += j * coefficients[j] * std::pow(z, j - 1);
    }
    for (int j = 2; j <= degree; ++j) {
        y_double_prime += j * (j - 1) * coefficients[j] * std::pow(z, j - 2);
    }

    // Now y_prime and y_double_prime contain the first and second derivatives at the point z

